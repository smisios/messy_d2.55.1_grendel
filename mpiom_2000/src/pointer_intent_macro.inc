#ifdef HAVE_POINTER_WITH_ATTR_INTENT
#  define POINTERINTENT(intent) POINTER,INTENT(intent)
#else
#  define POINTERINTENT(intent) POINTER
#endif

