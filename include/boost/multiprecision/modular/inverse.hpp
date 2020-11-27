//---------------------------------------------------------------------------//
// Copyright (c) 2020 Mikhail Komarov <nemo@nil.foundation>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef BOOST_MULTIPRECISION_MODULAR_BACKENDS_INVERSE_HPP
#define BOOST_MULTIPRECISION_MODULAR_BACKENDS_INVERSE_HPP

#include <boost/container/vector.hpp>

#include <boost/type_traits/is_integral.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_int/cpp_int_config.hpp>

namespace boost {
namespace multiprecision {
namespace backends {

template <typename Backend>
Backend eval_extended_euclidean_algorithm(Backend& a, Backend& b, Backend& x, Backend& y)
{
   if (eval_is_zero(a))
   {
      typedef typename mpl::front<typename Backend::unsigned_types>::type ui_type;
      x = ui_type(0u);
      y = ui_type(1u);
      return b;
   }
   Backend x1, y1, tmp = b;
   eval_modulus(tmp, a);
   Backend d = eval_extended_euclidean_algorithm(tmp, a, x1, y1);
   tmp       = b;
   eval_divide(tmp, a);
   eval_multiply(tmp, x1);
   x = y1;
   eval_subtract(x, tmp);
   y = x1;
   return d;
}

template <typename Backend>
Backend eval_inverse_extended_euclidean_algorithm(const Backend& a, const Backend& m)
{
   Backend aa = a, mm = m, x, y, g;
   typedef typename mpl::front<typename Backend::unsigned_types>::type ui_type;
   g = eval_extended_euclidean_algorithm(aa, mm, x, y);
   if (!eval_eq(g, ui_type(1u)))
   {
      return ui_type(0u);
   }
   else
   {
      eval_modulus(x, m);
      eval_add(x, m);
      eval_modulus(x, m);
      return x;
   }
}


template <typename Backend>
typename mpl::front<typename Backend::signed_types>::type eval_monty_inverse(typename mpl::front<typename Backend::signed_types>::type a)
{
   typedef typename mpl::front<typename Backend::signed_types>::type si_type;

   return eval_monty_inverse<si_type>(a);
}

template <typename T>
T eval_monty_inverse(T a)
{
   if (a % 2 == 0)
   {
      throw std::invalid_argument("monty_inverse only valid for odd integers");
   }

   T b = 1;
   T r = 0;

   for (size_t i = 0; i != sizeof(T) * CHAR_BIT; ++i)
   {
      const T bi = b % 2;
      r >>= 1;
      r += bi << (sizeof(T) * CHAR_BIT - 1);

      b -= a * bi;
      b >>= 1;
   }

   r = (~static_cast<T>(0) - r) + 1;

   return r;
}

template <typename Backend>
void eval_monty_inverse(Backend& res, const Backend& a, const Backend& p, const Backend& k)
{

   using default_ops::eval_modulus;
   using default_ops::eval_subtract;
   using default_ops::eval_abs;
   using default_ops::eval_gt;

   typedef typename mpl::front<typename Backend::unsigned_types>::type ui_type;
   Backend                                                             zero = ui_type(0u);
   Backend                                                             one  = ui_type(1u);
   Backend                                                             two  = ui_type(2u);
   Backend c, tmp;

   c = eval_inverse_extended_euclidean_algorithm(a, p);

   Backend bi = one, bt, i = zero, k_negone = k, xi, nextp = one;
   eval_subtract(k_negone, one);
   res = zero;

   ui_type kn = cpp_int(k_negone);

   while (!eval_eq(i, k))
   {
      xi = bi;
      eval_multiply(xi, c);
      eval_modulus(xi, p);

      if (eval_get_sign(xi) < 0)
      {
         tmp = xi;
         eval_abs(tmp, tmp);
         eval_modulus(tmp, p);
         xi = p;
         eval_subtract(xi, tmp);
      }

      tmp = a;
      eval_multiply(tmp, xi);
      eval_subtract(bi, tmp);
      eval_divide(bi, p);

      tmp = xi;
      eval_multiply(tmp, nextp);
      eval_multiply(nextp, p);
      eval_add(res, tmp);
      eval_add(i, one);
   }
}

template <typename Backend>
Backend eval_inverse_mod_pow2(Backend& a1, size_t k)
{
   typedef typename mpl::front<typename Backend::unsigned_types>::type ui_type;

   if (eval_integer_modulus(a1, 2) == 0)
      return 0;

   Backend b0 = 1, xi, bi, x = 0, a_xi;

   for (std::size_t i = 0; i != k; ++i)
   {
      bi = b0;
      eval_right_shift(bi, 1);
      xi = bi;
      eval_multiply(a_xi, a1, xi);
      eval_subtract(b0, b0, a_xi);
      eval_right_shift(b0, 1);
      eval_left_shift(x, 1);
      eval_add(x, xi);
   }

   return x;
}

}
}
} // namespace boost::multiprecision::backends

#endif

