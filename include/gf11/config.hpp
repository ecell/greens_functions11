#ifndef GF11_CONFIG_HPP
#define GF11_CONFIG_HPP

#if defined(GF11_HEADER_ONLY)
#  ifndef GF11_INLINE
#  define GF11_INLINE inline
#  endif
#else // pre-built
#  ifndef GF11_INLINE
#  define GF11_INLINE
#  endif
#endif

#endif//GF11_CONFIG_HPP
