 /* vim: set shiftwidth=2 softtabstop=2: */
/*
 *  ssss version 0.5  -  Copyright 2005,2006 B. Poettering
 * 
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 *  02111-1307 USA
 */

/*
 * http://point-at-infinity.org/ssss/
 *
 * This is an implementation of Shamir's Secret Sharing Scheme. See 
 * the project's homepage http://point-at-infinity.org/ssss/ for more 
 * information on this topic.
 *
 * This code links against the GNU multiprecision library "libgmp".
 * I compiled the code successfully with gmp 4.1.4.
 * You will need a system that has a /dev/random entropy source.
 *
 * Compile with 
 * "gcc -O2 -lgmp -o ssss-split ssss.c && ln ssss-split ssss-combine"
 *
 * Compile with -DNOMLOCK to obtain a version without memory locking.
 *
 * Report bugs to: ssss AT point-at-infinity.org
 *  
 */

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <termios.h>
#include <sys/mman.h>

#include <gmp.h>
#include <openssl/evp.h>

#ifdef __APPLE__
#include <CommonCrypto/CommonRandom.h>
#endif /* __APPLE__ */

#include "ssss.h"

#define DEBUG 0

/* coefficients of some irreducible polynomials over GF(2) */
static const unsigned char irred_coeff[] = {
  4,3,1,5,3,1,4,3,1,7,3,2,5,4,3,5,3,2,7,4,2,4,3,1,10,9,3,9,4,2,7,6,2,10,9,
  6,4,3,1,5,4,3,4,3,1,7,2,1,5,3,2,7,4,2,6,3,2,5,3,2,15,3,2,11,3,2,9,8,7,7,
  2,1,5,3,2,9,3,1,7,3,1,9,8,3,9,4,2,8,5,3,15,14,10,10,5,2,9,6,2,9,3,2,9,5,
  2,11,10,1,7,3,2,11,2,1,9,7,4,4,3,1,8,3,1,7,4,1,7,2,1,13,11,6,5,3,2,7,3,2,
  8,7,5,12,3,2,13,10,6,5,3,2,5,3,2,9,5,2,9,7,2,13,4,3,4,3,1,11,6,4,18,9,6,
  19,18,13,11,3,2,15,9,6,4,3,1,16,5,2,15,14,6,8,5,2,15,11,2,11,6,2,7,5,3,8,
  3,1,19,16,9,11,9,6,15,7,6,13,4,3,14,13,3,13,6,3,9,5,2,19,13,6,19,10,3,11,
  6,5,9,2,1,14,3,2,13,3,1,7,5,4,11,9,8,11,6,5,23,16,9,19,14,6,23,10,2,8,3,
  2,5,4,3,9,6,4,4,3,2,13,8,6,13,11,1,13,10,3,11,6,5,19,17,4,15,14,7,13,9,6,
  9,7,3,9,7,1,14,3,2,11,8,2,11,6,4,13,5,2,11,5,1,11,4,1,19,10,3,21,10,6,13,
  3,1,15,7,5,19,18,10,7,5,3,12,7,2,7,5,1,14,9,6,10,3,2,15,13,12,12,11,9,16,
  9,7,12,9,3,9,5,2,17,10,6,24,9,3,17,15,13,5,4,3,19,17,8,15,6,3,19,6,1 };

unsigned int degree;
mpz_t poly;
int cprng;
struct termios echo_orig, echo_off;

#define mpz_lshift(A, B, l) mpz_mul_2exp(A, B, l)
#define mpz_sizeinbits(A) (mpz_cmp_ui(A, 0) ? mpz_sizeinbase(A, 2) : 0)

void cprng_read_buf(void *buf, size_t buf_len);

/* emergency abort and warning functions */

void fatal(char *msg)
{
  tcsetattr(0, TCSANOW, &echo_orig);
  fprintf(stderr, "%sFATAL: %s.\n", isatty(2) ? "\a" : "", msg);
  exit(1);
}

void warning(opts_t opts, char *msg)
{
  if ((opts & opts_quiet) == 0) {
    fprintf(stderr, "%sWARNING: %s.\n", isatty(2) ? "\a" : "", msg);
  }
}

/* field arithmetic routines */

int field_size_valid(int deg)
{
  return (deg >= 8) && (deg <= MAXDEGREE) && (deg % 8 == 0);
}

/* initialize 'poly' to a bitfield representing the coefficients of an
   irreducible polynomial of degree 'deg' */

void field_init(int deg)
{
  assert(field_size_valid(deg));
  mpz_init_set_ui(poly, 0);
  mpz_setbit(poly, deg);
  mpz_setbit(poly, irred_coeff[3 * (deg / 8 - 1) + 0]);
  mpz_setbit(poly, irred_coeff[3 * (deg / 8 - 1) + 1]);
  mpz_setbit(poly, irred_coeff[3 * (deg / 8 - 1) + 2]);
  mpz_setbit(poly, 0);
  degree = deg;
}

void field_deinit(void)
{
  mpz_clear(poly);
}

/* I/O routines for GF(2^deg) field elements */

int field_import(mpz_t x, const char *s, size_t s_len, int hexmode)
{
  int result = -1;
  if (hexmode) {
    require_action(s_len <= degree / 4, out, warning(opts_none, "input string too long"));
    if (s_len < degree / 4) {
      warning(opts_none, "input string too short, adding null padding on the left");
    }
    require_action(!mpz_set_str(x, s, 16) && mpz_cmp_ui(x, 0) >= 0, out, warning(opts_none, "invalid syntax"));
  } else {
    require_action(s_len <= degree / 8, out, warning(opts_none, "input string too long"));
    mpz_import(x, s_len, 1, 1, 0, 0, s);
  }

  result = 0;
out:
  return result;
}

#if DEBUG
void dump_bytes(const char *label, const void *buf, size_t buf_len)
{
    size_t i;
    const unsigned char *b = (const unsigned char *)buf;

    fprintf(stderr, "%s(%zd):", label, buf_len);
    for (i = 0; i < buf_len; i++) {
        fprintf(stderr, "%02x", b[i]);
    }
    fprintf(stderr, "\n");
}
#endif

void prompt_for_input(const char *prompt, char *buf, size_t buf_len, bool echo)
{
  if (!echo) {
    tcsetattr(0, TCSANOW, &echo_off);
  }

  if (prompt) {
    fprintf(stderr, "%s", prompt);
  }

  if (!fgets(buf, buf_len, stdin)) {
    fatal("I/O error while reading prompted input");
  }

  /* replace newline with nul */
  buf[strcspn(buf, "\r\n")] = '\0';

  if (!echo) {
    fprintf(stderr, "\n");
  }

  if (!echo) {
    tcsetattr(0, TCSANOW, &echo_orig);
  }
}

typedef struct {
  unsigned char version;
  union {
    struct {
      uint32_t pbkdf_iterations;
      unsigned char pbkdf_salt[PKBDF_SALT_LEN];
      unsigned char aes_gcm_iv[GCM_IV_LEN];
      unsigned char aes_gcm_tag[GCM_TAG_LEN];
    } v0;
  };
} __attribute__((packed)) enc_params_t;

/* encrypt or decrypt the secret */
int encrypt_or_decrypt_secret(opts_t opts, bool encrypt, const char *buf, size_t buf_len, const char *passphrase, size_t passphrase_len, int pbkdf_iterations, enc_params_t *params, char *result_out)
{
  int result = -1;
  unsigned char enc_key[32] = {};
  unsigned char ciphertext[MAXLINELEN] = {};
  int ciphertext_len = (int)sizeof(ciphertext);
  unsigned char plaintext[MAXLINELEN] = {};
  int plaintext_len = (int)sizeof(plaintext);
  int len;
  EVP_CIPHER_CTX *ctx = NULL;

  assert(buf_len < INT_MAX);

  if (encrypt) {
    /* generate new encryption params: pbkdf salt and gcm iv */
    params->version = 0;
    params->v0.pbkdf_iterations = pbkdf_iterations;
    cprng_read_buf(params->v0.pbkdf_salt, sizeof(params->v0.pbkdf_salt));
    cprng_read_buf(params->v0.aes_gcm_iv, sizeof(params->v0.aes_gcm_iv));
  }

  /* derive AES key from encryption passphrase */
  switch(params->version) {
  case 0:
    if (1 != PKCS5_PBKDF2_HMAC(passphrase, (int)passphrase_len,
	params->v0.pbkdf_salt, (int)sizeof(params->v0.pbkdf_salt),
	params->v0.pbkdf_iterations, EVP_sha256(), (int)sizeof(enc_key), enc_key)) {
      fatal("pbkdf2 failed");
    }
#if DEBUG
    fprintf(stderr, "passphrase %.*s\n", (int)passphrase_len, passphrase);
    dump_bytes("enc_key", enc_key, sizeof(enc_key));
#endif

    if (NULL == (ctx = EVP_CIPHER_CTX_new())) {
      fatal("encryption context failed");
    }

    if (encrypt) {
      if (1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, enc_key, params->v0.aes_gcm_iv)) {
	fatal("encryption init failed");
      }

      if (1 != EVP_EncryptUpdate(ctx, ciphertext, &len, (unsigned char *)buf, buf_len)) {
	fatal("encryption update failed");
      }
      ciphertext_len = len;

      if (1 != EVP_EncryptFinal_ex(ctx, ciphertext + len, &len)) {
	fatal("encryption finalization failed");
      }
      ciphertext_len += len;

      if (1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_GET_TAG, GCM_TAG_LEN,
	  params->v0.aes_gcm_tag)) {
	fatal("fetching encryption tag failed");
      }

#if DEBUG
      fprintf(stderr, "plaintext len: %zu, ciphertext_len: %d\n", buf_len, ciphertext_len);
#endif
      if (ciphertext_len != (int)buf_len) {
	fatal("bug");
      }
      memcpy(result_out, ciphertext, ciphertext_len);
#if DEBUG
      dump_bytes("\nciphertext", result_out, ciphertext_len);
#endif
    } else {
      if (1 != EVP_DecryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, enc_key, params->v0.aes_gcm_iv)) {
	fatal("decryption init failed");
      }
      if (1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_TAG, GCM_TAG_LEN, params->v0.aes_gcm_tag)) {
	fatal("decryption tag failed");
      }
      if (1 != EVP_DecryptUpdate(ctx, plaintext, &len, (unsigned char *)buf, buf_len)) {
	fatal("decryption update failed");
      }
      plaintext_len = len;
      require_action(1 == EVP_DecryptFinal_ex(ctx, plaintext + len, &len), out, warning(opts, "decryption finalization failed -- wrong passcode?"));
      plaintext_len += len;
      if (plaintext_len != (int)buf_len) {
	fatal("bug");
      }
      memcpy(result_out, plaintext, plaintext_len);
    }
    break;
  default:
    fatal("invalid version");
  }

  result = 0;
out:
  if (ctx) {
    EVP_CIPHER_CTX_free(ctx);
  }
  /* clear out sensitive stuff from the stack */
  memset_s(enc_key, sizeof(enc_key), 0, sizeof(enc_key));
  memset_s(plaintext, sizeof(plaintext), 0, sizeof(plaintext));
  memset_s(ciphertext, sizeof(ciphertext), 0, sizeof(ciphertext));
  return result;
}

typedef int (*printer_f)(void *, const char * restrict, ...);

void _binary_print(printer_f p, void *printer_ctx, void *binary, size_t binary_len, bool hexmode)
{
  bool warn = false, printable = false;
  unsigned char *b = (unsigned char *)binary;
  for (unsigned int i = 0; i < binary_len; i++) {
      if (hexmode) {
          (*p)(printer_ctx, "%02x", b[i]);
      } else {
          printable = (b[i] >= 32) && (b[i] < 127);
          warn = warn || ! printable;
          (*p)(printer_ctx, "%c", printable ? b[i] : '.');
      }
  }
  if (warn && p == (printer_f)fprintf)
      warning(opts_none, "binary data detected, use -x mode instead");
}

void binary_print(FILE* stream, void *binary, size_t binary_len, bool hexmode)
{
  _binary_print((printer_f)fprintf, stream, binary, binary_len, hexmode);
}

void binary_append(fmt_buf_t *ctx, void *binary, size_t binary_len, bool hexmode)
{
  _binary_print((printer_f)fmt_buf_append, ctx, binary, binary_len, hexmode);
}

void enc_params_append(fmt_buf_t *ctx, bool hexmode, enc_params_t *params)
{
  /*  print version */
  binary_append(ctx, &params->version, sizeof(params->version), hexmode);
  /* print versioned stuct */
  switch(params->version) {
  case 0:
    binary_append(ctx, &params->v0, sizeof(params->v0), hexmode);
    break;
  default:
    break;
  }
}

void field_append(fmt_buf_t *ctx, const mpz_t x, bool hexmode)
{
  char buf[MAXDEGREE / 8 + 1] = {};
  size_t t;
  mpz_export(buf, &t, 1, 1, 0, 0, x);
  assert(t <= sizeof(buf));
  if (hexmode) {
    for (int i = degree / 8 - t; i; i--) {
      fmt_buf_append(ctx, "00");
    }
  }
  binary_append(ctx, buf, t, hexmode);
}

/* basic field arithmetic in GF(2^deg) */

void field_add(mpz_t z, const mpz_t x, const mpz_t y)
{
  mpz_xor(z, x, y);
}

void field_mult(mpz_t z, const mpz_t x, const mpz_t y)
{
  mpz_t b;
  unsigned int i;
  assert(z != y);
  mpz_init_set(b, x);
  if (mpz_tstbit(y, 0))
    mpz_set(z, b);
  else
    mpz_set_ui(z, 0);
  for(i = 1; i < degree; i++) {
    mpz_lshift(b, b, 1);
    if (mpz_tstbit(b, degree))
      mpz_xor(b, b, poly);
    if (mpz_tstbit(y, i))
      mpz_xor(z, z, b);
  }
  mpz_clear(b);
}

void field_invert(mpz_t z, const mpz_t x)
{
  mpz_t u, v, g, h;
  int i;
  assert(mpz_cmp_ui(x, 0));
  mpz_init_set(u, x);
  mpz_init_set(v, poly);
  mpz_init_set_ui(g, 0);
  mpz_set_ui(z, 1);
  mpz_init(h);
  while (mpz_cmp_ui(u, 1)) {
    i = mpz_sizeinbits(u) - mpz_sizeinbits(v);
    if (i < 0) {
      mpz_swap(u, v);
      mpz_swap(z, g);
      i = -i;
    }
    mpz_lshift(h, v, i);
    mpz_xor(u, u, h);
    mpz_lshift(h, g, i);
    mpz_xor(z, z, h);
  }
  mpz_clear(u); mpz_clear(v); mpz_clear(g); mpz_clear(h);
}

/* routines for the random number generator */

void cprng_init(void)
{
#ifdef __APPLE__
#else
  if ((cprng = open(RANDOM_SOURCE, O_RDONLY)) < 0)
    fatal("couldn't open " RANDOM_SOURCE);
#endif
}

void cprng_deinit(void)
{
#ifdef __APPLE__
#else
  if (close(cprng) < 0)
    fatal("couldn't close " RANDOM_SOURCE);
#endif
}

void cprng_read_buf(void *buf, size_t buf_len)
{
#ifdef __APPLE__
  if (kCCSuccess != CCRandomGenerateBytes(buf, buf_len)) {
    fatal("couldn't generated random bytes");
  }
#else
  if ((ssize_t)buf_len != read(cprng, buf, buf_len)) {
    fatal("couldn't read from " RANDOM_SOURCE);
  }
#endif
}

void cprng_read(mpz_t x)
{
  char buf[MAXDEGREE / 8];
#ifdef __APPLE__
  if (kCCSuccess != CCRandomGenerateBytes(buf, sizeof(buf))) {
    fatal("couldn't generated random bytes");
  }
#else
  unsigned int count;
  int i;
  for(count = 0; count < degree / 8; count += i)
    if ((i = read(cprng, buf + count, degree / 8 - count)) < 0) {
      close(cprng);
      fatal("couldn't read from " RANDOM_SOURCE);
    }
#endif
  mpz_import(x, degree / 8, 1, 1, 0, 0, buf);
}

/* a 64 bit pseudo random permutation (based on the XTEA cipher) */

void encipher_block(uint32_t *v) 
{
  uint32_t sum = 0, delta = 0x9E3779B9;
  int i;
  for(i = 0; i < 32; i++) {
    v[0] += (((v[1] << 4) ^ (v[1] >> 5)) + v[1]) ^ sum;
    sum += delta;
    v[1] += (((v[0] << 4) ^ (v[0] >> 5)) + v[0]) ^ sum;
  }
}

void decipher_block(uint32_t *v)
{
  uint32_t sum = 0xC6EF3720, delta = 0x9E3779B9;
  int i;
  for(i = 0; i < 32; i++) {
    v[1] -= ((v[0] << 4 ^ v[0] >> 5) + v[0]) ^ sum;
    sum -= delta;
    v[0] -= ((v[1] << 4 ^ v[1] >> 5) + v[1]) ^ sum;
  }
}

void encode_slice(uint8_t *data, int idx, int len, 
		  void (*process_block)(uint32_t*))
{
  uint32_t v[2];
  int i;
  for(i = 0; i < 2; i++)
    v[i] = data[(idx + 4 * i) % len] << 24 | 
      data[(idx + 4 * i + 1) % len] << 16 | 
      data[(idx + 4 * i + 2) % len] << 8 | data[(idx + 4 * i + 3) % len];
  process_block(v);
  for(i = 0; i < 2; i++) {
    data[(idx + 4 * i + 0) % len] = v[i] >> 24;
    data[(idx + 4 * i + 1) % len] = (v[i] >> 16) & 0xff;
    data[(idx + 4 * i + 2) % len] = (v[i] >> 8) & 0xff;
    data[(idx + 4 * i + 3) % len] = v[i] & 0xff;
  }
}

enum encdec {ENCODE, DECODE};

void encode_mpz(mpz_t x, enum encdec encdecmode)
{
  uint8_t v[(MAXDEGREE + 8) / 16 * 2];
  size_t t;
  int i;
  memset(v, 0, (degree + 8) / 16 * 2);
  mpz_export(v, &t, -1, 2, 1, 0, x);
  if (degree % 16 == 8)
    v[degree / 8 - 1] = v[degree / 8];
  if (encdecmode == ENCODE)             /* 40 rounds are more than enough!*/
    for(i = 0; i < 40 * ((int)degree / 8); i += 2)
      encode_slice(v, i, degree / 8, encipher_block);
  else
    for(i = 40 * (degree / 8) - 2; i >= 0; i -= 2)
      encode_slice(v, i, degree / 8, decipher_block);
  if (degree % 16 == 8) {
    v[degree / 8] = v[degree / 8 - 1];
    v[degree / 8 - 1] = 0;
  }
  mpz_import(x, (degree + 8) / 16, -1, 2, 1, 0, v);
  assert(mpz_sizeinbits(x) <= degree);
}

/* evaluate polynomials efficiently */

void horner(int n, mpz_t y, const mpz_t x, const mpz_t coeff[])
{
  int i;
  mpz_set(y, x);
  for(i = n - 1; i; i--) {
    field_add(y, y, coeff[i]);
    field_mult(y, y, x);
  }
  field_add(y, y, coeff[0]);
}

/* calculate the secret from a set of shares solving a linear equation system */

#define MPZ_SWAP(A, B) \
  do { mpz_set(h, A); mpz_set(A, B); mpz_set(B, h); } while(0)

int restore_secret(int n, mpz_t (*A)[n], mpz_t b[]) 
{
  mpz_t (*AA)[n] = (mpz_t (*)[n])A;
  int i, j, k, found;
  mpz_t h;
  mpz_init(h);
  for(i = 0; i < n; i++) {
    if (! mpz_cmp_ui(AA[i][i], 0)) {
      for(found = 0, j = i + 1; j < n; j++)
	if (mpz_cmp_ui(AA[i][j], 0)) {
	  found = 1;
	  break;
	}
      if (! found) 
	return -1;
      for(k = i; k < n; k++) 
	MPZ_SWAP(AA[k][i], AA[k][j]);
      MPZ_SWAP(b[i], b[j]);
    }
    for(j = i + 1; j < n; j++) {
      if (mpz_cmp_ui(AA[i][j], 0)) {
	for(k = i + 1; k < n; k++) {
	  field_mult(h, AA[k][i], AA[i][j]);
	  field_mult(AA[k][j], AA[k][j], AA[i][i]);
	  field_add(AA[k][j], AA[k][j], h);
	}
	field_mult(h, b[i], AA[i][j]);
	field_mult(b[j], b[j], AA[i][i]);
	field_add(b[j], b[j], h);
      }
    }
  }
  field_invert(h, AA[n - 1][n - 1]);
  field_mult(b[n - 1], b[n - 1], h);
  mpz_clear(h);
  return 0;
}

/* split a secret into shares, optionally encrypted using passphrase */
int split(opts_t opts, int threshold, int numbner, int security, const char *token, char *secret, int secret_len, const char *passphrase, size_t passphrase_len, int pbkdf_iterations, char **shares_out[])
{
  int result = -1;
  char **shares = NULL;
  unsigned int fmt_len;
  mpz_t x, y, coeff[threshold];
  int i;
  bool do_encrypt = passphrase != NULL;
  bool do_token = token != NULL;
  enc_params_t params = {};
  fmt_buf_t ctx;
  char ciphertext[secret_len];

  for(fmt_len = 1, i = numbner; i >= 10; i /= 10, fmt_len++);

  cprng_init();
  if (do_encrypt) {
    require_noerr_action(encrypt_or_decrypt_secret(opts, true, secret, secret_len, passphrase, passphrase_len,
	  pbkdf_iterations, &params, ciphertext), out, warning(opts, "encryption failed"));
    secret = ciphertext;
  }

  if (! security) {
    security = (opts & opts_hex) ? 4 * ((secret_len + 1) & ~1): 8 * secret_len;
  }

  require_action(field_size_valid(security), out, warning(opts, "security level invalid (secret too long?)"));

  if ((opts & opts_quiet) != 0)
      fprintf(stderr, "Using a %d bit security level.\n", security);

  field_init(security);

  mpz_init(coeff[0]);
  require_noerr_action(field_import(coeff[0], secret, secret_len, (opts & opts_hex) != 0), out, warning(opts, "import failed"));

  if ((opts & opts_diffusion) != 0) {
    if (degree >= 64)
      encode_mpz(coeff[0], ENCODE);
    else
      warning(opts, "security level too small for the diffusion layer");
  }

  for(i = 1; i < threshold; i++) {
    mpz_init(coeff[i]);
    cprng_read(coeff[i]);
  }
  cprng_deinit();

  mpz_init(x);
  mpz_init(y);
  require_action(shares = malloc(numbner * sizeof(char*)), out, warning(opts, "malloc"));
  for (i = 0; i < numbner; i++) {
    mpz_set_ui(x, i + 1);
    horner(threshold, y, x, (const mpz_t*)coeff);
    fmt_buf_init(&ctx);
    if (do_token || do_encrypt) {
      if (do_token) {
	fmt_buf_append(&ctx, "%s", token);
      }
      if (do_encrypt) {
	fmt_buf_append(&ctx, ":");
	enc_params_append(&ctx, true, &params);
      }
      fmt_buf_append(&ctx, "-");
    }
    fmt_buf_append(&ctx, "%0*d-", fmt_len, i + 1);
    field_append(&ctx, y, true);
    shares[i] = fmt_buf_get(&ctx);
  }
  mpz_clear(x);
  mpz_clear(y);

  for(i = 0; i < threshold; i++)
    mpz_clear(coeff[i]);

  *shares_out = shares;

  result = 0;

out:
  fmt_buf_destroy(&ctx);
  field_deinit();

  return result;
}

#define HEX2BIN(ch) \
  (((ch) >= '0' && (ch) <= '9') ? (ch) - '0' : \
  ((ch) >= 'A' && (ch) <= 'F') ? (ch) - 'A' + 10 : \
  ((ch) >= 'a' && (ch) <= 'f') ? (ch) - 'a' + 10 : 0)

void decode_enc_params(const char *encoded_params, enc_params_t *params)
{
  size_t len = strlen(encoded_params);
  if ((len & 1) == 1) {
    warning(opts_none, "invalid enc params for share, ignoring");
  } else {
    /* decode the first byte, to get the vesion */
    params->version = (HEX2BIN(encoded_params[0])) << 4 | (HEX2BIN(encoded_params[1]));
    encoded_params += 2;
    switch (params->version) {
    case 0:
      if (len/2 - sizeof(params->version) != sizeof(params->v0)) {
	fatal("invalid length of enc params");
      }
      for (int i = 0; i < (int)sizeof(params->v0); i++) {
	((char *)&params->v0)[i] = (HEX2BIN(encoded_params[i*2])) << 4 | (HEX2BIN(encoded_params[i*2+1]));
      }
      break;
    default:
      fatal("invalid encryption params version");
      break;
    }
  }
#if DEBUG
  dump_bytes("params", params, sizeof(*params));
#endif
}

/* Prompt for shares, calculate the secret */

int combine(opts_t opts, char shares[][MAXLINELEN], int threshold, const char *passphrase, size_t passphrase_len, char **secret)
{
  int result = -1;
  mpz_t A[threshold][threshold], y[threshold], x;
  char *a, *b, *c = NULL;
  int i, j;
  unsigned s = 0;
  enc_params_t params = {};
  char buf[MAXDEGREE / 8 + 1] = {};
  char plaintext[MAXDEGREE / 8 + 1] = {};
  fmt_buf_t ctx;
  size_t t;

  assert(threshold >= 2);

  mpz_init(x);
  for (i = 0; i < threshold; i++) {
    require_action(a = strchr(shares[i], '-'), out, warning(opts, "invalid syntax"));
    *a++ = 0;
    if ((b = strchr(a, '-'))) {
      *b++ = 0;
      c = shares[i];
    } else {
      b = a, a = shares[i];
    }

#if DEBUG
    fprintf(stderr, "a -> %s\n", a);
    fprintf(stderr, "b -> %s\n", b);
    fprintf(stderr, "c -> %s\n", c);
#endif
    /* if our token has enc params, decode them */
    if (c && (c = strchr(c, ':'))) {
      c++;
      decode_enc_params(c, &params);
    }

    if (! s) {
      s = 4 * strlen(b);
      require_action(field_size_valid(s), out, warning(opts, "share has illegal length"));
      field_init(s);
    }
    else
      require_action(s == 4 * strlen(b), out, warning(opts, "shares have different security levels"));

    require_action(j = atoi(a), out, warning(opts, "invalid share"));
    mpz_set_ui(x, j);
    mpz_init_set_ui(A[threshold - 1][i], 1);
    for(j = threshold - 2; j >= 0; j--) {
      mpz_init(A[j][i]);
      field_mult(A[j][i], A[j + 1][i], x);
    }
    mpz_init(y[i]);
    require_noerr_action(field_import(y[i], b, strlen(b), 1), out, warning(opts, "field import"));
    field_mult(x, x, A[0][i]);
    field_add(y[i], y[i], x);
  }
  mpz_clear(x);
  require_noerr_action(restore_secret(threshold, A, y), out, warning(opts, "shares inconsistent. Perhaps a single share was used twice"));

  if ((opts & opts_diffusion) != 0) {
    if (degree >= 64)
      encode_mpz(y[threshold - 1], DECODE);
    else
      if ((opts & opts_quiet) == 0)
	warning(opts, "security level too small for the diffusion layer");
  }

  mpz_export(buf, &t, 1, 1, 0, 0, y[threshold - 1]);

  if ((opts & opts_encrypt) != 0) {
    require_noerr_action(encrypt_or_decrypt_secret(opts, false, buf, t, passphrase, passphrase_len, 0, &params, plaintext),
	out, warning(opts, "encryption failed"));
    memcpy(buf, plaintext, t);
  }
  fmt_buf_init(&ctx);
  binary_append(&ctx, buf, t, (opts & opts_hex) != 0);
  *secret = fmt_buf_get(&ctx);

  for (i = 0; i < threshold; i++) {
    for (j = 0; j < threshold; j++)
      mpz_clear(A[i][j]);
    mpz_clear(y[i]);
  }

  result = 0;
out:
  field_deinit();
  return result;
}

void combine_shares(opts_t opts, int threshold)
{
  char buf[threshold][MAXLINELEN];
  char encryption_passphrase[128];
  char *secret = NULL;
  int i;

  if ((opts & opts_quiet) == 0)
    printf("Enter %d shares separated by newlines:\n", threshold);

  for (i = 0; i < threshold; i++) {
    if ((opts & opts_quiet) == 0)
      printf("Share [%d/%d]: ", i + 1, threshold);

    prompt_for_input(NULL, buf[i], sizeof(buf[i]), true);
  }

  if ((opts & opts_encrypt) != 0) {
    /* prompt for encryption passphrase */
    prompt_for_input("Enter the encryption passphrase that protects the secret: ",
	encryption_passphrase, sizeof(encryption_passphrase), false);
  }

  combine(opts, buf, threshold,
      (opts & opts_encrypt) != 0 ? encryption_passphrase : NULL,
      (opts & opts_encrypt) != 0 ? strlen(encryption_passphrase) : 0,
      &secret);

  if ((opts & opts_quiet) == 0)
    fprintf(stderr, "Resulting secret: ");
  fprintf(stderr, "%s\n", secret);
  free(secret);
  memset_s(encryption_passphrase, sizeof(encryption_passphrase), 0, sizeof(encryption_passphrase));
}

void dump_fmt_buf_state(fmt_buf_t *ctx)
{
  printf("::\n");
  printf("capacity: %ld\n", ctx->end - ctx->start);
  printf("used: %ld\n", ctx->next_unused - ctx->start);
  printf("::n\n");
}

void fmt_buf_init(fmt_buf_t *ctx)
{
  ctx->start = ctx->next_unused = ctx->end = NULL;
}


void fmt_buf_destroy(fmt_buf_t *ctx)
{
  if (ctx->start) {
    free(ctx->start);
  }
  fmt_buf_init(ctx);
}

char *fmt_buf_get(fmt_buf_t *ctx)
{
  if (ctx->start) {
    return strdup(ctx->start);
  } else {
    return NULL;
  }
}

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define INIT_ALLOC_SIZE 64
void fmt_buf_append(fmt_buf_t *ctx, const char * restrict format, ...)
{
  va_list ap;
  if (ctx->start == NULL) {
    /* uninitialized */
    if (NULL == (ctx->next_unused = ctx->start = calloc(1, INIT_ALLOC_SIZE))) {
      fatal("malloc");
    }
    ctx->end = ctx->start + INIT_ALLOC_SIZE;
  }

  va_start(ap, format);
  int n = vsnprintf(ctx->next_unused, ctx->end - ctx->next_unused, format, ap);
  va_end(ap);
  if (n < 0) {
    fatal("vsnpritnf");
  }
  if (n >= (ctx->end - ctx->next_unused)) {
    /* truncated, resize */
    size_t old_offset = ctx->next_unused - ctx->start;
    size_t new_size =  MAX(2*(ctx->end - ctx->start), (ctx->end - ctx->start) + n);
    if (NULL == (ctx->start = realloc(ctx->start, new_size))) {
      fatal("malloc");
    }
    ctx->next_unused = ctx->start + old_offset;
    ctx->end = ctx->start + new_size;

    va_start(ap, format);
    n = vsnprintf(ctx->next_unused, ctx->end - ctx->next_unused, format, ap);
    if (n < 0 || n >= (ctx->end - ctx->next_unused)) {
      fatal("bug");
    }
    va_end(ap);
  }
  ctx->next_unused += n;
  assert(ctx->next_unused <= ctx->end);
  assert(*ctx->next_unused == '\0');
}


#ifndef TEST
int main(int argc, char *argv[])
{
  char *name;
  int i;
  opts_t opts = opts_diffusion;
  int opt_security = 0;
  int opt_threshold = -1;
  int opt_number = -1;
  int opt_pbkdf_iterations = PBKDF_ITERATIONS;
  char *opt_token = NULL;

#if ! NOMLOCK
  if (mlockall(MCL_CURRENT | MCL_FUTURE) < 0)
    switch(errno) {
    case ENOMEM:
      warning(opts_none, "couldn't get memory lock (ENOMEM, try to adjust RLIMIT_MEMLOCK!)");
      break;
    case EPERM:
      warning(opts_none, "couldn't get memory lock (EPERM, try UID 0!)");
      break;
    case ENOSYS:
      warning(opts_none, "couldn't get memory lock (ENOSYS, kernel doesn't allow page locking)");
      break;
    default:
      warning(opts_none, "couldn't get memory lock");
      break;
    }
#endif

  if (getuid() != geteuid())
    seteuid(getuid());

  tcgetattr(0, &echo_orig);
  echo_off = echo_orig;
  echo_off.c_lflag &= ~ECHO;

  opts |= (argc == 1) ? opts_help : 0;
  while((i = getopt(argc, argv, "vDEI:hqQxs:t:n:w:")) != -1)
    switch(i) {
    case 'v': opts |= opts_showversion; break;
    case 'h': opts |= opts_help; break;
    case 'q': opts |= opts_quiet; break;
    case 'Q': opts |= opts_QUIET|opts_quiet; break;
    case 'x': opts |= opts_hex; break;
    case 's': opt_security = atoi(optarg); break;
    case 't': opt_threshold = atoi(optarg); break;
    case 'n': opt_number = atoi(optarg); break;
    case 'w': opt_token = optarg; break;
    case 'D': opts &= ~opts_diffusion; break;
    case 'E': opts |= opts_encrypt; opts &= ~opts_diffusion; break;
    case 'I': opt_pbkdf_iterations = atoi(optarg); break;
    default:
      exit(1);
    }
  if ((opts & opts_help) == 0 && (argc != optind))
    fatal("invalid argument");

  if ((name = strrchr(argv[0], '/')) == NULL)
    name = argv[0];

  if (strstr(name, "split")) {
    if ((opts & (opts_help|opts_showversion))) {
      puts("Split secrets using Shamir's Secret Sharing Scheme.\n"
	   "\n"
	   "ssss-split -t threshold -n shares [-w token] [-s level]"
	   " [-x] [-q] [-Q] [-D] [-E] [-I iterations] [-v]"
	   );
      if (opts & opts_showversion)
	puts("\nVersion: " VERSION);
      exit(0);
    }

    if (opt_threshold < 2)
      fatal("invalid parameters: invalid threshold value");

    if (opt_number < opt_threshold)
      fatal("invalid parameters: number of shares smaller than threshold");

    if (opt_security && ! field_size_valid(opt_security))
      fatal("invalid parameters: invalid security level");

    if (opt_token && (strlen(opt_token) > MAXTOKENLEN))
      fatal("invalid parameters: token too long");

    int deg;
    size_t secret_len;
    char buf[MAXLINELEN];
    char encryption_passphrase[128] = {};
    if ((opts & opts_quiet) == 0) {
      printf("Generating shares using a (%d,%d) scheme with ", 
	  opt_threshold, opt_number);
      if (opt_security)
	printf("a %d bit", opt_security);
      else
	printf("dynamic");
      printf(" security level.\n");

      deg = opt_security ? opt_security : MAXDEGREE;
      fprintf(stderr, "Enter the secret, ");
      if ((opts & opts_hex) != 0)
	fprintf(stderr, "as most %d hex digits: ", deg / 4);
      else
	fprintf(stderr, "at most %d ASCII characters: ", deg / 8);
    }
    prompt_for_input(NULL, buf, sizeof(buf), false);
    secret_len = strlen(buf);

    if (opts & opts_encrypt) {
      /* prompt for encryption passphrase */
      prompt_for_input("Enter the encryption passphrase that protects the secret: ",
	  encryption_passphrase, sizeof(encryption_passphrase), false);
    }

    char **shares = NULL;

    split(opts, opt_threshold, opt_number, opt_security, opt_token, buf, (int)secret_len,
	(opts & opts_encrypt) ? encryption_passphrase : NULL,
	(opts & opts_encrypt) ? strlen(encryption_passphrase) : 0,
	(opts & opts_encrypt) ? opt_pbkdf_iterations : 0, &shares);

    for (int i = 0; i < opt_number; i ++){
      fprintf(stdout, "%s\n", shares[i]);
      free(shares[i]);
    }
    free(shares);
    memset_s(encryption_passphrase, sizeof(encryption_passphrase), 0, sizeof(encryption_passphrase));
  }
  else {
    if ((opts & (opts_help|opts_showversion))) {
      puts("Combine shares using Shamir's Secret Sharing Scheme.\n"
	   "\n"
	   "ssss-combine -t threshold [-x] [-q] [-Q] [-D] [-E] [-v]");
      if (opts & opts_showversion)
	puts("\nVersion: " VERSION);
      exit(0);
    }

    if (opt_threshold < 2)
      fatal("invalid parameters: invalid threshold value");

    combine_shares(opts, opt_threshold);
  }
  return 0;
}
#endif /* TEST */
