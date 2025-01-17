//---------------------------------------------------------------------------//
// Copyright (c) 2020 Mikhail Komarov <nemo@nil.foundation>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef BOOST_MULTIPRECISION_INVERSE_HPP
#define BOOST_MULTIPRECISION_INVERSE_HPP

#include <boost/container/vector.hpp>

#include <boost/type_traits/is_integral.hpp>

#include <nil/crypto3/multiprecision/cpp_int.hpp>
#include <nil/crypto3/multiprecision/cpp_int/cpp_int_config.hpp>
#include <nil/crypto3/multiprecision/modular/modular_adaptor.hpp>
#include <nil/crypto3/multiprecision/modular/modular_adaptor_fixed.hpp>
#include <nil/crypto3/multiprecision/modular/inverse.hpp>

namespace nil {
    namespace crypto3 {
        namespace multiprecision {

            template<typename Backend, expression_template_option ExpressionTemplates>
            constexpr number<Backend, ExpressionTemplates>
                inverse_extended_euclidean_algorithm(const number<Backend, ExpressionTemplates>& n,
                                                     const number<Backend, ExpressionTemplates>& mod) {
                number<Backend, ExpressionTemplates> result;
                backends::eval_inverse_extended_euclidean_algorithm(result.backend(), n.backend(), mod.backend());
                return result;
            }

            template<typename Backend, expression_template_option ExpressionTemplates>
            constexpr number<modular_adaptor<Backend>, ExpressionTemplates> inverse_extended_euclidean_algorithm(
                const number<modular_adaptor<Backend>, ExpressionTemplates>& modular) {
                number<Backend, ExpressionTemplates> new_base, res;
                number<modular_adaptor<Backend>, ExpressionTemplates> res_mod;

                modular.backend().mod_data().adjust_regular(new_base.backend(), modular.backend().base_data());
                backends::eval_inverse_extended_euclidean_algorithm(res.backend(), new_base.backend(), modular.backend().mod_data().get_mod().backend());
                assign_components(res_mod.backend(), res.backend(), modular.backend().mod_data().get_mod().backend());

                return res_mod;
            }

            template<typename Backend, expression_template_option ExpressionTemplates>
            constexpr number<Backend, ExpressionTemplates>
                monty_inverse(const number<Backend, ExpressionTemplates>& a,
                              const number<Backend, ExpressionTemplates>& p,
                              const number<Backend, ExpressionTemplates>& k) {
                number<Backend, ExpressionTemplates> res;
                backends::eval_monty_inverse(res.backend(), a.backend(), p.backend(), k.backend());
                return res;
            }

            template<typename Backend, expression_template_option ExpressionTemplates>
            constexpr number<Backend, ExpressionTemplates>
            inverse_mod(const number<Backend, ExpressionTemplates>& a,
                          const number<Backend, ExpressionTemplates>& p) {
                number<Backend, ExpressionTemplates> res;
                backends::eval_inverse_mod(res.backend(), a.backend(), p.backend());
                return res;
            }

            template<typename Backend, expression_template_option ExpressionTemplates>
            constexpr number<modular_adaptor<Backend>, ExpressionTemplates>
            inverse_mod(const number<modular_adaptor<Backend>, ExpressionTemplates>& modular) {
                number<modular_adaptor<Backend>, ExpressionTemplates> res;
                backends::eval_inverse_mod(res.backend(), modular.backend());
                return res;
            }

        }    // namespace multiprecision
    }        // namespace crypto3
}    // namespace nil

#endif
