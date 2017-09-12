/*
# This file is part of the Astrometry.net suite.
# Licensed under a 3-clause BSD style license - see LICENSE
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

#include <assert.h>

#include "log.h"

static log_t g_logger;

static log_t* get_logger() {
    return &g_logger;
}

void log_init_structure(log_t* logger, enum log_level level) {
	logger->level = level;
    logger->f = stdout;
	logger->timestamp = FALSE;
	logger->t0 = time();
	logger->logfunc = NULL;
	logger->baton = NULL;
}

void log_init(enum log_level level) {
	log_init_structure(get_logger(), level);
}

void log_set_level(enum log_level level) {
    get_logger()->level = level;
}

void log_set_timestamp(anbool b) {
	get_logger()->timestamp = b;
}

void log_to(FILE* fid) {
	get_logger()->f = fid;
}

void log_to_fd(int fd) {
    // MEMLEAK
    FILE* fid = fdopen(fd, "a");
    log_to(fid);
}

void log_use_function(logfunc_t func, void* baton) {
	log_t* l = get_logger();
	l->logfunc = func;
	l->baton = baton;
}

log_t* log_create(enum log_level level) {
	log_t* logger = calloc(1, sizeof(log_t));
	return logger;
}

void log_free(log_t* log) {
	assert(log);
	free(log);
}

static void loglvl(const log_t* logger, enum log_level level,
				   const char* file, int line, const char* func,
                   const char* format, va_list va) {
	if (level > logger->level)
		return;
	if (logger->f) {
		if (logger->timestamp)
			//fprintf(logger->f, "[%6i: %.3f] ", (int)getpid(), time() - logger->t0);
		//fprintf(logger->f, "%s:%i ", file, line);
		vfprintf(logger->f, format, va);
		fflush(logger->f);
	}
	if (logger->logfunc) {
		logger->logfunc(logger->baton, level, file, line, func, format, va);
	}
}

void log_loglevel(enum log_level level,
				  const char* file, int line, const char* func,
				  const char* format, ...) {
    va_list va;
    va_start(va, format);
    loglvl(get_logger(), level, file, line, func, format, va);
    va_end(va);
}

int log_get_level() {
    return get_logger()->level;
}

FILE* log_get_fid() {
	return get_logger()->f;
}

#define LOGGER_TEMPLATE(name, level)									\
	void																\
	name##_(const log_t* logger, const char* file, int line, const char* func, const char* format, ...) { \
		va_list va;														\
		va_start(va, format);											\
		loglvl(logger, level, file, line, func, format, va);				\
		va_end(va);														\
	}																	\
	void																\
	name(const char* file, int line, const char* func, const char* format, ...) { \
		va_list va;														\
		va_start(va, format);											\
		loglvl(get_logger(), level, file, line, func, format, va);		\
		va_end(va);														\
	}																	\

LOGGER_TEMPLATE(log_logerr,  LOG_ERROR);
LOGGER_TEMPLATE(log_logmsg,  LOG_MSG);
LOGGER_TEMPLATE(log_logverb, LOG_VERB);
LOGGER_TEMPLATE(log_logdebug,LOG_ALL);
