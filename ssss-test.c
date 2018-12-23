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
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <assert.h>

#include "ssss.h"

int fmt_buf_test(void)
{
  fmt_buf_t ctx;
  fmt_buf_init(&ctx);

  //dump_fmt_buf_state(&ctx);
  fmt_buf_append(&ctx, "0123456789abcde%c0123456789abcde%c", 'f', 'f');
  assert(ctx.next_unused - ctx.start == 32);
  //dump_fmt_buf_state(&ctx);
  //printf("%s\n", ctx.start);
  fmt_buf_append(&ctx, "0123456789abcde%c0123456789abcd%c", 'f', 'e');
  assert(ctx.next_unused - ctx.start == 63);
  //dump_fmt_buf_state(&ctx);
  //printf("%s\n", ctx.start);

  fmt_buf_append(&ctx, "%c", 'f');
  assert(ctx.next_unused - ctx.start == 64);
  //dump_fmt_buf_state(&ctx);
  //printf("%s\n", ctx.start);

  fmt_buf_append(&ctx, "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcde%c", 'f');
  assert(ctx.next_unused - ctx.start == 208);
  //dump_fmt_buf_state(&ctx);
  //printf("%s\n", ctx.start);

  for (int i = 0; i < 1000; i++) {
    fmt_buf_append(&ctx, "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcde%03d", i);
  }
  assert(ctx.next_unused - ctx.start == 208 + (1000*146));
  //dump_fmt_buf_state(&ctx);
  //printf("%s\n", ctx.start);
  fmt_buf_destroy(&ctx);

  return 0;
}

int split_combine_test(opts_t opts, int secret_len, int token_len, int threshold)
{
  int result = -1;
  char in_shares[threshold][MAXLINELEN];
  char secret[secret_len + 1];
  char token[token_len + 1];
  int t1, t2;
  char **out_shares = NULL;
  char *result_out = NULL;
  int number = 10;

  memset(secret, 'A', secret_len);
  secret[secret_len] = '\0';
  if (token_len > 0) {
    memset(token, 'T', token_len);
    token[token_len] = '\0';
  }

  require_noerr(split(opts, threshold, number, 0, token_len > 0 ? token : NULL, secret, secret_len,
	(opts & opts_encrypt) ? "password" : NULL,
	(opts & opts_encrypt) ? strlen("password") : 0,
	(opts & opts_encrypt) ? 1 /* low iteration count for speed */ : 0,
	&out_shares), out);

  for (int i = 0; i < 50; i++) {
    /* pick two random and distinct shares */
    do {
      t1 = arc4random() % number;
      t2 = arc4random() % number;
    } while (t1 == t2);

    strlcpy(in_shares[0], out_shares[t1], sizeof(in_shares[0]));
    strlcpy(in_shares[1], out_shares[t2], sizeof(in_shares[1]));

    /* make sure they're combinable */
    if (opts & opts_encrypt) {
      /* try wrong passcode, assert combining fails */
      require(combine(opts|opts_quiet, in_shares, threshold,
	    (opts & opts_encrypt) ? "passworD" : NULL,
	    (opts & opts_encrypt) ? strlen("passworD") : 0,
	    &result_out), out);

      strlcpy(in_shares[0], out_shares[t1], sizeof(in_shares[0]));
      strlcpy(in_shares[1], out_shares[t2], sizeof(in_shares[1]));
    }

    /* try corrupting a share */
    size_t len = strlen(in_shares[0]);
    switch(in_shares[0][len - 1]) {
      case '0':
	in_shares[0][len - 1] = '1';
	break;
      default:
	in_shares[0][len - 1] = '0';
	break;
    }
    if (0 == combine(opts|opts_quiet, in_shares, threshold,
	  (opts & opts_encrypt) ? "password" : NULL,
	  (opts & opts_encrypt) ? strlen("password") : 0,
	  &result_out)) {
      /* assert secret is different */
      require(strcmp(result_out, secret) != 0, out);
      free(result_out);
      result_out = NULL;
    }

    strlcpy(in_shares[0], out_shares[t1], sizeof(in_shares[0]));
    strlcpy(in_shares[1], out_shares[t2], sizeof(in_shares[1]));

    require_noerr(combine(opts, in_shares, threshold,
	(opts & opts_encrypt) ? "password" : NULL,
	(opts & opts_encrypt) ? strlen("password") : 0,
	&result_out), out);
    if (strcmp(result_out, secret) != 0) {
      printf("FAILED: opts: %x, token: %s, iteration: %d, %d:%d\n", opts, token, i, t1, t2);
      printf("shares: %s, %s\n", in_shares[0], in_shares[1]);
      abort();
    }
    free(result_out);
    result_out = NULL;
  }

  /* free shares */
  for (int i = 0; i < threshold; i ++){
    free(out_shares[i]);
  }
  free(out_shares);

  result = 0;
out:
  return result;
}

void backward_compat_test(void)
{
  const char *answer = "ssss-0.1";
  char *result = NULL;
  int threshold = 2;
  char in_shares[3][MAXLINELEN] = {
    "1-188586b9f78418272ad22629f3ac654f",
    "2-310b0d73ef08304ec031d9c69008b8cb",
    "3-298e8bca188c286999908c9c4e94f3b5",
  };
  opts_t opts = 0;

  if (combine(opts, in_shares, threshold, NULL, 0, &result)) {
    fprintf(stderr, "test failed");
  }
  assert(strcmp(result, answer) == 0);
  free(result);
}

int main(int argc, char *argv[])
{
  (void)argc;
  (void)argv;

  struct {
    opts_t opts;
    int threshold;
  } tests[] = {
    { opts_quiet|opts_none,      2 },
    { opts_quiet|opts_none,      2 },
    { opts_quiet|opts_diffusion, 2 },
    { opts_quiet|opts_diffusion, 2 },
    { opts_quiet|opts_encrypt,   2 },
    { opts_quiet|opts_encrypt,   2 },
  };

  int secret_lengths[] = {
    1, 16, 128,
  };

  int token_lengths[] = {
    0, 1, 16,
  };

#define SIZ(x) sizeof(x)/sizeof(x[0])

  fmt_buf_test();
  backward_compat_test();
  for (size_t i = 0; i < SIZ(tests); i++) {
    for (size_t s = 0; s < SIZ(secret_lengths); s++) {
      for (size_t t = 0; t < SIZ(token_lengths); t++) {
	printf("testing tokenlen: %d, secretlen %d, opts 0x%x, threshold: %d", token_lengths[t], secret_lengths[s], tests[i].opts, tests[i].threshold);
	int r = split_combine_test(tests[i].opts, secret_lengths[s], token_lengths[t], tests[i].threshold);
	printf(" result = %s\e[0m\n", r == 0 ? "\e[1;32mOK" : "\e[1;31mFAIL");
      }
    }
  }

  printf("ok\n");
  return EXIT_SUCCESS;
}
