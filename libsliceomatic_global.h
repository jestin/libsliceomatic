#ifndef LIBSLICEOMATIC_GLOBAL_H
#define LIBSLICEOMATIC_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(LIBSLICEOMATIC_LIBRARY)
#  define LIBSLICEOMATICSHARED_EXPORT Q_DECL_EXPORT
#else
#  define LIBSLICEOMATICSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // LIBSLICEOMATIC_GLOBAL_H
