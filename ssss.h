#ifndef __SSSS_H
#define __SSSS_H

#define VERSION "0.5"
#define RANDOM_SOURCE "/dev/random"
#define MAXDEGREE 1024
#define MAXTOKENLEN 128
#define MAXLINELEN (MAXTOKENLEN + 1 + 10 + 1 + MAXDEGREE / 4 + 10)
#define PBKDF_ITERATIONS 100000
#define PKBDF_SALT_LEN 8
#define GCM_IV_LEN 16
#define GCM_TAG_LEN 16

typedef struct {
  char *start;
  char *next_unused;
  char *end;
} fmt_buf_t;

typedef enum {
  opts_none = 0,
  opts_hex = 1 << 0,
  opts_diffusion = 1 << 1,
  opts_quiet = 1 << 2,
  opts_QUIET = 1 << 3,
  opts_encrypt = 1 << 4,
  opts_help = 1 << 5,
  opts_showversion = 1 << 6,
} opts_t;

void fmt_buf_init(fmt_buf_t *ctx);
void fmt_buf_destroy(fmt_buf_t *ctx);
void fmt_buf_append(fmt_buf_t *ctx, const char * restrict format, ...);
char *fmt_buf_get(fmt_buf_t *ctx);
void dump_fmt_buf_state(fmt_buf_t *ctx);
void fatal(char *msg);

int split(opts_t opts, int threshold, int numbner, int security, const char *token, char *secret, int secret_len, const char *passphrase, size_t passphrase_len, int pbkdf_iterations, char **shares_out[]);
int combine(opts_t opts, char shares[][MAXLINELEN], int threshold, const char *passphrase, size_t passphrase_len, char **secret);

#define likely(x)   __builtin_expect((x), 1)
#define unlikely(x) __builtin_expect((x), 0)

#define require(thing, label) do { if (unlikely(!(thing))) { goto label; } } while (0);
#define require_noerr(thing, label) do { if (unlikely((thing))) { goto label; } } while (0);
#define require_action(thing, label, action) do { if (unlikely(!(thing))) { {action;} goto label; } } while (0);
#define require_noerr_action(thing, label, action) do { if (unlikely((thing))) { {action;} goto label; } } while (0);

#endif /* __SSSS_H */
