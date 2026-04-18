// See LICENSE for licensing info.
// deflate.c -- RFC 1951 raw DEFLATE compress + decompress.
//
// Adapted from cute_png.h's cp_inflate / cp_lz77_compress + associated
// helpers (https://github.com/RandyGaul/cute_headers). Same author as this
// codebase; MIT-licensed in original form.
//
// Differences from cute_png.h:
//   - Stripped all PNG-specific machinery (IHDR/IDAT/IEND chunks, filtering,
//     CRC, Adler32, zlib wrapper). Produces + consumes raw DEFLATE streams.
//   - Prefixed helpers with dfl_ instead of cp_ to make intent clear inside
//     nudge.
//   - Output buffer uses CK_ALLOC / CK_REALLOC so allocator hooks propagate.
//
// Public API:
//   int deflate_compress  (const void* in, int in_len, uint8_t** out, int* out_len);
//   int deflate_decompress(const void* in, int in_len, void* out, int out_len);
//
// deflate_compress emits a single fixed-Huffman block (BFINAL=1, BTYPE=01).
// Acceptable for save files where one-shot compression is fine -- no need
// for zlib wrapping, no need to stream.

// ---------------------------------------------------------------------------
// DEFLATE tables (RFC 1951 §3.2.5 / §3.2.6 / §3.2.7)

#define DFL_LOOKUP_BITS   9
#define DFL_LOOKUP_COUNT  (1 << DFL_LOOKUP_BITS)

static uint8_t  dfl_fixed_table[288 + 32] = {
	8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
	8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
	8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
	7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
};
static uint8_t  dfl_permutation_order[19] = { 16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15 };
static uint8_t  dfl_len_extra_bits[29 + 2] = { 0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0,  0,0 };
static uint32_t dfl_len_base[29 + 2] = { 3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258,  0,0 };
static uint8_t  dfl_dist_extra_bits[30 + 2] = { 0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,  0,0 };
static uint32_t dfl_dist_base[30 + 2] = { 1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577, 0,0 };

// ---------------------------------------------------------------------------
// Inflate (decompress).

typedef struct dfl_in_t
{
	uint64_t bits;
	int count;
	uint32_t* words;
	int word_count;
	int word_index;
	int bits_left;

	int final_word_available;
	uint32_t final_word;

	char* out;
	char* out_end;
	char* begin;

	uint16_t lookup[DFL_LOOKUP_COUNT];
	uint32_t lit[288];
	uint32_t dst[32];
	uint32_t len[19];
	uint32_t nlit;
	uint32_t ndst;
	uint32_t nlen;
} dfl_in_t;

static int dfl_would_overflow(dfl_in_t* s, int num_bits)
{
	return (s->bits_left + s->count) - num_bits < 0;
}

static char* dfl_ptr(dfl_in_t* s)
{
	assert(!(s->bits_left & 7));
	return (char*)(s->words + s->word_index) - (s->count / 8);
}

static uint64_t dfl_peak_bits(dfl_in_t* s, int num_bits_to_read)
{
	if (s->count < num_bits_to_read) {
		if (s->word_index < s->word_count) {
			uint32_t word = s->words[s->word_index++];
			s->bits |= (uint64_t)word << s->count;
			s->count += 32;
		} else if (s->final_word_available) {
			uint32_t word = s->final_word;
			s->bits |= (uint64_t)word << s->count;
			s->count += s->bits_left;
			s->final_word_available = 0;
		}
	}
	return s->bits;
}

static uint32_t dfl_consume_bits(dfl_in_t* s, int num_bits_to_read)
{
	assert(s->count >= num_bits_to_read);
	uint32_t bits = (uint32_t)(s->bits & (((uint64_t)1 << num_bits_to_read) - 1));
	s->bits >>= num_bits_to_read;
	s->count -= num_bits_to_read;
	s->bits_left -= num_bits_to_read;
	return bits;
}

static uint32_t dfl_read_bits(dfl_in_t* s, int num_bits_to_read)
{
	assert(num_bits_to_read <= 32);
	assert(num_bits_to_read >= 0);
	assert(s->bits_left > 0);
	assert(s->count <= 64);
	assert(!dfl_would_overflow(s, num_bits_to_read));
	dfl_peak_bits(s, num_bits_to_read);
	return dfl_consume_bits(s, num_bits_to_read);
}

static uint32_t dfl_rev16(uint32_t a)
{
	a = ((a & 0xAAAA) >>  1) | ((a & 0x5555) << 1);
	a = ((a & 0xCCCC) >>  2) | ((a & 0x3333) << 2);
	a = ((a & 0xF0F0) >>  4) | ((a & 0x0F0F) << 4);
	a = ((a & 0xFF00) >>  8) | ((a & 0x00FF) << 8);
	return a;
}

static int dfl_build(dfl_in_t* s, uint32_t* tree, uint8_t* lens, int sym_count)
{
	int counts[16] = { 0 };
	int codes[16] = { 0 };
	int first[16] = { 0 };

	for (int i = 0; i < sym_count; ++i) ++counts[lens[i]];
	for (int n = 1; n < 16; ++n) {
		codes[n] = (codes[n - 1] + counts[n - 1]) << 1;
		first[n] = first[n - 1] + counts[n - 1];
	}

	if (s) memset(s->lookup, 0, sizeof(s->lookup));
	for (int i = 0; i < sym_count; ++i) {
		int len = lens[i];
		if (len != 0) {
			assert(len < 16);
			uint32_t code = codes[len]++;
			uint32_t slot = first[len]++;
			tree[slot] = (code << (32 - len)) | (i << 4) | len;

			if (s && len <= DFL_LOOKUP_BITS) {
				int j = dfl_rev16(code) >> (16 - len);
				while (j < (1 << DFL_LOOKUP_BITS)) {
					s->lookup[j] = (uint16_t)((len << DFL_LOOKUP_BITS) | i);
					j += (1 << len);
				}
			}
		}
	}
	return first[15];
}

static int dfl_stored(dfl_in_t* s)
{
	dfl_read_bits(s, s->count & 7);  // skip partial byte
	uint16_t LEN = (uint16_t)dfl_read_bits(s, 16);
	uint16_t NLEN = (uint16_t)dfl_read_bits(s, 16);
	if (LEN != (uint16_t)(~NLEN)) return 0;
	if (s->bits_left / 8 < (int)LEN) return 0;
	char* p = dfl_ptr(s);
	memcpy(s->out, p, LEN);
	s->out += LEN;
	return 1;
}

static int dfl_fixed(dfl_in_t* s)
{
	s->nlit = dfl_build(s, s->lit, dfl_fixed_table, 288);
	s->ndst = dfl_build(0, s->dst, dfl_fixed_table + 288, 32);
	return 1;
}

static int dfl_decode(dfl_in_t* s, uint32_t* tree, int hi)
{
	uint64_t bits = dfl_peak_bits(s, 16);
	uint32_t search = (dfl_rev16((uint32_t)bits) << 16) | 0xFFFF;
	int lo = 0;
	while (lo < hi) {
		int guess = (lo + hi) >> 1;
		if (search < tree[guess]) hi = guess;
		else lo = guess + 1;
	}
	uint32_t key = tree[lo - 1];
	uint32_t len = (32 - (key & 0xF));
	assert((search >> len) == (key >> len));
	(void)len;
	dfl_consume_bits(s, key & 0xF);
	return (key >> 4) & 0xFFF;
}

static int dfl_dynamic(dfl_in_t* s)
{
	uint8_t lenlens[19] = { 0 };
	int nlit = 257 + dfl_read_bits(s, 5);
	int ndst = 1 + dfl_read_bits(s, 5);
	int nlen = 4 + dfl_read_bits(s, 4);
	for (int i = 0 ; i < nlen; ++i)
		lenlens[dfl_permutation_order[i]] = (uint8_t)dfl_read_bits(s, 3);

	s->nlen = dfl_build(0, s->len, lenlens, 19);
	uint8_t lens[288 + 32];
	for (int n = 0; n < nlit + ndst;) {
		int sym = dfl_decode(s, s->len, s->nlen);
		switch (sym) {
		case 16: for (int i =  3 + dfl_read_bits(s, 2); i; --i, ++n) lens[n] = lens[n - 1]; break;
		case 17: for (int i =  3 + dfl_read_bits(s, 3); i; --i, ++n) lens[n] = 0; break;
		case 18: for (int i = 11 + dfl_read_bits(s, 7); i; --i, ++n) lens[n] = 0; break;
		default: lens[n++] = (uint8_t)sym; break;
		}
	}
	s->nlit = dfl_build(s, s->lit, lens, nlit);
	s->ndst = dfl_build(0, s->dst, lens + nlit, ndst);
	return 1;
}

static int dfl_block(dfl_in_t* s)
{
	while (1) {
		int symbol = dfl_decode(s, s->lit, s->nlit);
		if (symbol < 256) {
			if (s->out + 1 > s->out_end) return 0;
			*s->out++ = (char)symbol;
		} else if (symbol > 256) {
			symbol -= 257;
			int length = dfl_read_bits(s, dfl_len_extra_bits[symbol]) + dfl_len_base[symbol];
			int dsym = dfl_decode(s, s->dst, s->ndst);
			int back = dfl_read_bits(s, dfl_dist_extra_bits[dsym]) + dfl_dist_base[dsym];
			if (s->out - back < s->begin) return 0;
			if (s->out + length > s->out_end) return 0;
			char* src = s->out - back;
			char* dst = s->out;
			s->out += length;
			if (back == 1) { memset(dst, *src, length); }
			else while (length--) *dst++ = *src++;
		} else break;
	}
	return 1;
}

int deflate_decompress(const void* in, int in_len, void* out, int out_len)
{
	dfl_in_t* s = (dfl_in_t*)CK_ALLOC(sizeof(dfl_in_t));
	memset(s, 0, sizeof(*s));
	s->bits_left = in_len * 8;

	int first_bytes = (int)((((size_t)in + 3) & ~3) - (size_t)in);
	s->words = (uint32_t*)((char*)in + first_bytes);
	s->word_count = (in_len - first_bytes) / 4;
	int last_bytes = ((in_len - first_bytes) & 3);

	for (int i = 0; i < first_bytes; ++i)
		s->bits |= (uint64_t)(((const uint8_t*)in)[i]) << (i * 8);

	s->final_word_available = last_bytes ? 1 : 0;
	s->final_word = 0;
	for (int i = 0; i < last_bytes; i++)
		s->final_word |= ((const uint8_t*)in)[in_len - last_bytes + i] << (i * 8);

	s->count = first_bytes * 8;
	s->out = (char*)out;
	s->out_end = s->out + out_len;
	s->begin = (char*)out;

	int ok = 1;
	int bfinal;
	do {
		bfinal = dfl_read_bits(s, 1);
		int btype = dfl_read_bits(s, 2);
		switch (btype) {
		case 0: ok = dfl_stored(s); break;
		case 1: dfl_fixed(s); ok = dfl_block(s); break;
		case 2: dfl_dynamic(s); ok = dfl_block(s); break;
		case 3: ok = 0; break;
		}
		if (!ok) break;
	} while (!bfinal);

	CK_FREE(s);
	return ok;
}

// ---------------------------------------------------------------------------
// Deflate (compress). Fixed-Huffman single-block encoder with LZ77.

typedef struct dfl_lz77_t
{
	int head[4096];         // hash table heads (-1 = empty)
	int prev[32768];        // chain links (32K sliding window)
	const uint8_t* window;
	int window_size;
} dfl_lz77_t;

typedef struct dfl_out_t
{
	uint8_t* data;
	int len;
	int cap;
	uint32_t bits;  // bit accumulator, sentinel-marked (0x80 = empty)
} dfl_out_t;

static void dfl_put8(dfl_out_t* s, uint32_t a)
{
	if (s->len >= s->cap) {
		s->cap = s->cap ? s->cap * 2 : 256;
		s->data = (uint8_t*)CK_REALLOC(s->data, (size_t)s->cap);
	}
	s->data[s->len++] = (uint8_t)a;
}

static void dfl_put_bits(dfl_out_t* s, uint32_t data, uint32_t bitcount)
{
	while (bitcount--) {
		uint32_t prev = s->bits;
		s->bits = (s->bits >> 1) | ((data & 1) << 7);
		data >>= 1;
		if (prev & 1) {
			dfl_put8(s, s->bits);
			s->bits = 0x80;
		}
	}
}

// put bits MSB first
static void dfl_put_bitsr(dfl_out_t* s, uint32_t data, uint32_t bitcount)
{
	while (bitcount--)
		dfl_put_bits(s, data >> bitcount, 1);
}

static uint32_t dfl_lz77_hash(const uint8_t* p)
{
	return ((uint32_t)(p[0] << 8) ^ (uint32_t)(p[1] << 4) ^ (uint32_t)p[2]) & 0xFFF;
}

static void dfl_lz77_init(dfl_lz77_t* lz, const uint8_t* data, int size)
{
	memset(lz->head, 0xFF, sizeof(lz->head));
	memset(lz->prev, 0xFF, sizeof(lz->prev));
	lz->window = data;
	lz->window_size = size;
}

static void dfl_lz77_insert(dfl_lz77_t* lz, int pos)
{
	if (pos + 2 >= lz->window_size) return;
	uint32_t hash = dfl_lz77_hash(lz->window + pos);
	int slot = pos & 0x7FFF;
	lz->prev[slot] = lz->head[hash];
	lz->head[hash] = pos;
}

static int dfl_lz77_find_match(dfl_lz77_t* lz, int pos, int max_len, int* out_distance)
{
	const uint8_t* data = lz->window;
	int best_len = 2;
	int best_dist = 0;
	int chain_count = 0;
	const int max_chain = 256;

	if (max_len < 3) return 0;
	if (max_len > 258) max_len = 258;

	uint32_t hash = dfl_lz77_hash(data + pos);
	int match_pos = lz->head[hash];
	if (match_pos < 0) return 0;

	while (match_pos >= 0 && chain_count < max_chain) {
		int dist = pos - match_pos;
		if (dist > 32768 || dist <= 0) break;
		if (data[match_pos] == data[pos] && data[match_pos + best_len] == data[pos + best_len]) {
			int len = 0;
			while (len < max_len && data[match_pos + len] == data[pos + len]) len++;
			if (len > best_len) {
				best_len = len;
				best_dist = dist;
				if (len >= max_len) break;
			}
		}
		match_pos = lz->prev[match_pos & 0x7FFF];
		chain_count++;
	}
	if (best_len >= 3) { *out_distance = best_dist; return best_len; }
	return 0;
}

static void dfl_encode_literal(dfl_out_t* s, uint32_t v)
{
	     if (v < 144) dfl_put_bitsr(s, 0x030 + v -   0, 8);
	else if (v < 256) dfl_put_bitsr(s, 0x190 + v - 144, 9);
	else if (v < 280) dfl_put_bitsr(s, 0x000 + v - 256, 7);
	else              dfl_put_bitsr(s, 0x0c0 + v - 280, 8);
}

static void dfl_encode_distance(dfl_out_t* s, int distance)
{
	int code = 0;
	for (int i = 29; i >= 0; i--) {
		if (distance >= (int)dfl_dist_base[i]) { code = i; break; }
	}
	dfl_put_bitsr(s, code, 5);
	if (dfl_dist_extra_bits[code] > 0)
		dfl_put_bits(s, distance - dfl_dist_base[code], dfl_dist_extra_bits[code]);
}

static void dfl_encode_match(dfl_out_t* s, int length, int distance)
{
	int code = 0;
	for (int i = 28; i >= 0; i--) {
		if (length >= (int)dfl_len_base[i]) { code = i; break; }
	}
	dfl_encode_literal(s, 257 + code);
	if (dfl_len_extra_bits[code] > 0)
		dfl_put_bits(s, length - dfl_len_base[code], dfl_len_extra_bits[code]);
	dfl_encode_distance(s, distance);
}

static void dfl_lz77_compress(dfl_out_t* s, dfl_lz77_t* lz, const uint8_t* data, int len)
{
	dfl_lz77_init(lz, data, len);
	int pos = 0;
	while (pos < len) {
		int distance;
		int match_len = dfl_lz77_find_match(lz, pos, len - pos, &distance);
		if (match_len >= 3) {
			dfl_encode_match(s, match_len, distance);
			for (int i = 0; i < match_len; i++) dfl_lz77_insert(lz, pos + i);
			pos += match_len;
		} else {
			dfl_encode_literal(s, data[pos]);
			dfl_lz77_insert(lz, pos);
			pos++;
		}
	}
}

int deflate_compress(const void* in, int in_len, uint8_t** out, int* out_len)
{
	dfl_out_t s = { 0 };
	s.bits = 0x80;
	s.cap = 256;
	s.data = (uint8_t*)CK_ALLOC((size_t)s.cap);

	// Single fixed-Huffman block: BFINAL=1, BTYPE=01. 3 bits total, LSB-first.
	dfl_put_bits(&s, 3, 3);

	// LZ77 state is ~144 KB; heap-allocate to avoid blowing stacks.
	dfl_lz77_t* lz = (dfl_lz77_t*)CK_ALLOC(sizeof(dfl_lz77_t));
	dfl_lz77_compress(&s, lz, (const uint8_t*)in, in_len);
	CK_FREE(lz);

	// End-of-block symbol.
	dfl_encode_literal(&s, 256);
	// Pad to byte boundary (flush the bit accumulator).
	while (s.bits != 0x80) dfl_put_bits(&s, 0, 1);

	*out = s.data;
	*out_len = s.len;
	return 1;
}
