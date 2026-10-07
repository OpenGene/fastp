/* LD_PRELOAD shim for CI: report FAKE_NPROCS CPUs so fastp doesn't cap -w to the
   runner's core count (options.cpp), letting CI exercise high worker counts. */
#define _GNU_SOURCE
#include <stdlib.h>
#include <unistd.h>
#include <dlfcn.h>

static int fake(void) { const char* v = getenv("FAKE_NPROCS"); return v ? atoi(v) : 64; }
int get_nprocs(void) { return fake(); }
int get_nprocs_conf(void) { return fake(); }
long sysconf(int name) {
    static long (*real)(int) = 0;
    if (name == _SC_NPROCESSORS_ONLN || name == _SC_NPROCESSORS_CONF) return fake();
    if (!real) real = (long (*)(int))dlsym(RTLD_NEXT, "sysconf");
    return real(name);
}
