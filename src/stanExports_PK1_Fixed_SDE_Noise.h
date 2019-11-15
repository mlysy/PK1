// Generated by rstantools.  Do not edit by hand.

/*
    PK1 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PK1 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PK1.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.18.0

#include <stan/model/model_header.hpp>

namespace model_PK1_Fixed_SDE_Noise_namespace {

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;

static int current_statement_begin__;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_PK1_Fixed_SDE_Noise");
    reader.add_event(20, 20, "include", "/include/PK1_DE_Functions.stan");
    reader.add_event(20, 0, "start", "/include/PK1_DE_Functions.stan");
    reader.add_event(48, 28, "end", "/include/PK1_DE_Functions.stan");
    reader.add_event(48, 21, "restart", "model_PK1_Fixed_SDE_Noise");
    reader.add_event(104, 75, "end", "model_PK1_Fixed_SDE_Noise");
    return reader;
}

template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__>::type>::type, Eigen::Dynamic,1>
PK1_ODE(const T0__& D,
            const T1__& Cl,
            const T2__& Ka,
            const T3__& Ke,
            const std::vector<T4__>& t, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__>::type>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 23;
        local_scalar_t__ R;
        (void) R;  // dummy to suppress unused var warning

        stan::math::initialize(R, DUMMY_VAR__);
        stan::math::fill(R,DUMMY_VAR__);


        current_statement_begin__ = 24;
        stan::math::assign(R, ((((Ke * Ka) * D) / Cl) / (Ke - Ka)));
        current_statement_begin__ = 25;
        return stan::math::promote_scalar<fun_return_scalar_t__>(multiply(R,subtract(stan::math::exp(multiply(-(Ka),to_vector(t))),stan::math::exp(multiply(-(Ke),to_vector(t))))));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct PK1_ODE_functor__ {
    template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__>::type>::type, Eigen::Dynamic,1>
    operator()(const T0__& D,
            const T1__& Cl,
            const T2__& Ka,
            const T3__& Ke,
            const std::vector<T4__>& t, std::ostream* pstream__) const {
        return PK1_ODE(D, Cl, Ka, Ke, t, pstream__);
    }
};

template <bool propto, typename T0__, typename T1__, typename T2__, typename T3__, typename T4__, typename T5__, typename T6__, typename T7__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__, T7__>::type>::type
PK1_SDE_lpdf(const std::vector<T0__>& Xt,
                 const T1__& Ka,
                 const T2__& Ke,
                 const T3__& Cl,
                 const T4__& sigmaP,
                 const T5__& Dose,
                 const std::vector<T6__>& t,
                 const std::vector<T7__>& dt, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__, T7__>::type>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 31;
        local_scalar_t__ R;
        (void) R;  // dummy to suppress unused var warning

        stan::math::initialize(R, DUMMY_VAR__);
        stan::math::fill(R,DUMMY_VAR__);
        current_statement_begin__ = 32;
        local_scalar_t__ ll;
        (void) ll;  // dummy to suppress unused var warning

        stan::math::initialize(ll, DUMMY_VAR__);
        stan::math::fill(ll,DUMMY_VAR__);
        current_statement_begin__ = 33;
        local_scalar_t__ rho;
        (void) rho;  // dummy to suppress unused var warning

        stan::math::initialize(rho, DUMMY_VAR__);
        stan::math::fill(rho,DUMMY_VAR__);
        current_statement_begin__ = 34;
        local_scalar_t__ lambda;
        (void) lambda;  // dummy to suppress unused var warning

        stan::math::initialize(lambda, DUMMY_VAR__);
        stan::math::fill(lambda,DUMMY_VAR__);
        current_statement_begin__ = 35;
        local_scalar_t__ tau;
        (void) tau;  // dummy to suppress unused var warning

        stan::math::initialize(tau, DUMMY_VAR__);
        stan::math::fill(tau,DUMMY_VAR__);
        current_statement_begin__ = 36;
        int nObs(0);
        (void) nObs;  // dummy to suppress unused var warning

        stan::math::fill(nObs, std::numeric_limits<int>::min());


        current_statement_begin__ = 37;
        stan::math::assign(nObs, num_elements(t));
        current_statement_begin__ = 38;
        stan::math::assign(ll, 0.0);
        current_statement_begin__ = 39;
        stan::math::assign(R, ((((Ke * Ka) * Dose) / Cl) / (Ke - Ka)));
        current_statement_begin__ = 40;
        for (int ii = 1; ii <= (nObs - 1); ++ii) {

            current_statement_begin__ = 41;
            stan::math::assign(rho, stan::math::exp((-(Ke) * get_base1(dt,ii,"dt",1))));
            current_statement_begin__ = 42;
            stan::math::assign(lambda, ((R * stan::math::exp((-(Ka) * get_base1(t,ii,"t",1)))) * (stan::math::exp((-(Ka) * get_base1(dt,ii,"dt",1))) - rho)));
            current_statement_begin__ = 43;
            stan::math::assign(tau, (sigmaP * stan::math::sqrt(((1 - pow(rho,2)) / (2.0 * Ke)))));
            current_statement_begin__ = 44;
            stan::math::assign(ll, (ll + normal_log(get_base1(Xt,(ii + 1),"Xt",1),((rho * get_base1(Xt,ii,"Xt",1)) + lambda),tau)));
        }
        current_statement_begin__ = 46;
        return stan::math::promote_scalar<fun_return_scalar_t__>(ll);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__, typename T5__, typename T6__, typename T7__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__, T7__>::type>::type
PK1_SDE_lpdf(const std::vector<T0__>& Xt,
                 const T1__& Ka,
                 const T2__& Ke,
                 const T3__& Cl,
                 const T4__& sigmaP,
                 const T5__& Dose,
                 const std::vector<T6__>& t,
                 const std::vector<T7__>& dt, std::ostream* pstream__) {
    return PK1_SDE_lpdf<false>(Xt,Ka,Ke,Cl,sigmaP,Dose,t,dt, pstream__);
}


struct PK1_SDE_lpdf_functor__ {
    template <bool propto, typename T0__, typename T1__, typename T2__, typename T3__, typename T4__, typename T5__, typename T6__, typename T7__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__, T7__>::type>::type
    operator()(const std::vector<T0__>& Xt,
                 const T1__& Ka,
                 const T2__& Ke,
                 const T3__& Cl,
                 const T4__& sigmaP,
                 const T5__& Dose,
                 const std::vector<T6__>& t,
                 const std::vector<T7__>& dt, std::ostream* pstream__) const {
        return PK1_SDE_lpdf(Xt, Ka, Ke, Cl, sigmaP, Dose, t, dt, pstream__);
    }
};

#include <stan_meta_header.hpp>
 class model_PK1_Fixed_SDE_Noise : public prob_grad {
private:
    int nObs;
    int nSub;
    vector<vector<double> > Yt;
    vector<vector<double> > t;
    vector<double> D;
    double sdDef;
    vector<vector<double> > dt;
public:
    model_PK1_Fixed_SDE_Noise(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }

    model_PK1_Fixed_SDE_Noise(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }

    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;

        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning

        current_statement_begin__ = -1;

        static const char* function__ = "model_PK1_Fixed_SDE_Noise_namespace::model_PK1_Fixed_SDE_Noise";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        try {
            current_statement_begin__ = 51;
            context__.validate_dims("data initialization", "nObs", "int", context__.to_vec());
            nObs = int(0);
            vals_i__ = context__.vals_i("nObs");
            pos__ = 0;
            nObs = vals_i__[pos__++];
            current_statement_begin__ = 52;
            context__.validate_dims("data initialization", "nSub", "int", context__.to_vec());
            nSub = int(0);
            vals_i__ = context__.vals_i("nSub");
            pos__ = 0;
            nSub = vals_i__[pos__++];
            current_statement_begin__ = 54;
            validate_non_negative_index("Yt", "nSub", nSub);
            validate_non_negative_index("Yt", "nObs", nObs);
            context__.validate_dims("data initialization", "Yt", "double", context__.to_vec(nSub,nObs));
            validate_non_negative_index("Yt", "nSub", nSub);
            validate_non_negative_index("Yt", "nObs", nObs);
            Yt = std::vector<std::vector<double> >(nSub,std::vector<double>(nObs,double(0)));
            vals_r__ = context__.vals_r("Yt");
            pos__ = 0;
            size_t Yt_limit_1__ = nObs;
            for (size_t i_1__ = 0; i_1__ < Yt_limit_1__; ++i_1__) {
                size_t Yt_limit_0__ = nSub;
                for (size_t i_0__ = 0; i_0__ < Yt_limit_0__; ++i_0__) {
                    Yt[i_0__][i_1__] = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 55;
            validate_non_negative_index("t", "nSub", nSub);
            validate_non_negative_index("t", "nObs", nObs);
            context__.validate_dims("data initialization", "t", "double", context__.to_vec(nSub,nObs));
            validate_non_negative_index("t", "nSub", nSub);
            validate_non_negative_index("t", "nObs", nObs);
            t = std::vector<std::vector<double> >(nSub,std::vector<double>(nObs,double(0)));
            vals_r__ = context__.vals_r("t");
            pos__ = 0;
            size_t t_limit_1__ = nObs;
            for (size_t i_1__ = 0; i_1__ < t_limit_1__; ++i_1__) {
                size_t t_limit_0__ = nSub;
                for (size_t i_0__ = 0; i_0__ < t_limit_0__; ++i_0__) {
                    t[i_0__][i_1__] = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 56;
            validate_non_negative_index("D", "nSub", nSub);
            context__.validate_dims("data initialization", "D", "double", context__.to_vec(nSub));
            validate_non_negative_index("D", "nSub", nSub);
            D = std::vector<double>(nSub,double(0));
            vals_r__ = context__.vals_r("D");
            pos__ = 0;
            size_t D_limit_0__ = nSub;
            for (size_t i_0__ = 0; i_0__ < D_limit_0__; ++i_0__) {
                D[i_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 58;
            context__.validate_dims("data initialization", "sdDef", "double", context__.to_vec());
            sdDef = double(0);
            vals_r__ = context__.vals_r("sdDef");
            pos__ = 0;
            sdDef = vals_r__[pos__++];

            // validate, data variables
            current_statement_begin__ = 51;
            check_greater_or_equal(function__,"nObs",nObs,1);
            current_statement_begin__ = 52;
            check_greater_or_equal(function__,"nSub",nSub,1);
            current_statement_begin__ = 54;
            current_statement_begin__ = 55;
            current_statement_begin__ = 56;
            for (int k0__ = 0; k0__ < nSub; ++k0__) {
                check_greater_or_equal(function__,"D[k0__]",D[k0__],0);
            }
            current_statement_begin__ = 58;
            check_greater_or_equal(function__,"sdDef",sdDef,0);
            // initialize data variables
            current_statement_begin__ = 62;
            validate_non_negative_index("dt", "nSub", nSub);
            validate_non_negative_index("dt", "(nObs - 1)", (nObs - 1));
            dt = std::vector<std::vector<double> >(nSub,std::vector<double>((nObs - 1),double(0)));
            stan::math::fill(dt,DUMMY_VAR__);

            current_statement_begin__ = 64;
            for (int jj = 1; jj <= (nObs - 1); ++jj) {

                current_statement_begin__ = 65;
                for (int ii = 1; ii <= nSub; ++ii) {

                    current_statement_begin__ = 66;
                    stan::model::assign(dt, 
                                stan::model::cons_list(stan::model::index_uni(ii), stan::model::cons_list(stan::model::index_uni(jj), stan::model::nil_index_list())), 
                                (get_base1(get_base1(t,ii,"t",1),(jj + 1),"t",2) - get_base1(get_base1(t,ii,"t",1),jj,"t",2)), 
                                "assigning variable dt");
                }
            }

            // validate transformed data
            current_statement_begin__ = 62;

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 73;
            ++num_params_r__;
            current_statement_begin__ = 74;
            ++num_params_r__;
            current_statement_begin__ = 75;
            ++num_params_r__;
            current_statement_begin__ = 78;
            ++num_params_r__;
            current_statement_begin__ = 79;
            ++num_params_r__;
            current_statement_begin__ = 82;
            validate_non_negative_index("Xt", "nSub", nSub);
            validate_non_negative_index("Xt", "nObs", nObs);
            num_params_r__ += nSub * nObs;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~model_PK1_Fixed_SDE_Noise() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        if (!(context__.contains_r("Cl")))
            throw std::runtime_error("variable Cl missing");
        vals_r__ = context__.vals_r("Cl");
        pos__ = 0U;
        context__.validate_dims("initialization", "Cl", "double", context__.to_vec());
        double Cl(0);
        Cl = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,Cl);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable Cl: ") + e.what());
        }

        if (!(context__.contains_r("Ka")))
            throw std::runtime_error("variable Ka missing");
        vals_r__ = context__.vals_r("Ka");
        pos__ = 0U;
        context__.validate_dims("initialization", "Ka", "double", context__.to_vec());
        double Ka(0);
        Ka = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,Ka);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable Ka: ") + e.what());
        }

        if (!(context__.contains_r("Ke")))
            throw std::runtime_error("variable Ke missing");
        vals_r__ = context__.vals_r("Ke");
        pos__ = 0U;
        context__.validate_dims("initialization", "Ke", "double", context__.to_vec());
        double Ke(0);
        Ke = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,Ke);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable Ke: ") + e.what());
        }

        if (!(context__.contains_r("sigmaP")))
            throw std::runtime_error("variable sigmaP missing");
        vals_r__ = context__.vals_r("sigmaP");
        pos__ = 0U;
        context__.validate_dims("initialization", "sigmaP", "double", context__.to_vec());
        double sigmaP(0);
        sigmaP = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,sigmaP);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable sigmaP: ") + e.what());
        }

        if (!(context__.contains_r("sigmaM")))
            throw std::runtime_error("variable sigmaM missing");
        vals_r__ = context__.vals_r("sigmaM");
        pos__ = 0U;
        context__.validate_dims("initialization", "sigmaM", "double", context__.to_vec());
        double sigmaM(0);
        sigmaM = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,sigmaM);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable sigmaM: ") + e.what());
        }

        if (!(context__.contains_r("Xt")))
            throw std::runtime_error("variable Xt missing");
        vals_r__ = context__.vals_r("Xt");
        pos__ = 0U;
        validate_non_negative_index("Xt", "nSub", nSub);
        validate_non_negative_index("Xt", "nObs", nObs);
        context__.validate_dims("initialization", "Xt", "double", context__.to_vec(nSub,nObs));
        std::vector<std::vector<double> > Xt(nSub,std::vector<double>(nObs,double(0)));
        for (int i1__ = 0U; i1__ < nObs; ++i1__)
            for (int i0__ = 0U; i0__ < nSub; ++i0__)
                Xt[i0__][i1__] = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < nSub; ++i0__)
            for (int i1__ = 0U; i1__ < nObs; ++i1__)
                try {
            writer__.scalar_unconstrain(Xt[i0__][i1__]);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable Xt: ") + e.what());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        typedef T__ local_scalar_t__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        try {
            // model parameters
            stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);

            local_scalar_t__ Cl;
            (void) Cl;  // dummy to suppress unused var warning
            if (jacobian__)
                Cl = in__.scalar_lb_constrain(0,lp__);
            else
                Cl = in__.scalar_lb_constrain(0);

            local_scalar_t__ Ka;
            (void) Ka;  // dummy to suppress unused var warning
            if (jacobian__)
                Ka = in__.scalar_lb_constrain(0,lp__);
            else
                Ka = in__.scalar_lb_constrain(0);

            local_scalar_t__ Ke;
            (void) Ke;  // dummy to suppress unused var warning
            if (jacobian__)
                Ke = in__.scalar_lb_constrain(0,lp__);
            else
                Ke = in__.scalar_lb_constrain(0);

            local_scalar_t__ sigmaP;
            (void) sigmaP;  // dummy to suppress unused var warning
            if (jacobian__)
                sigmaP = in__.scalar_lb_constrain(0,lp__);
            else
                sigmaP = in__.scalar_lb_constrain(0);

            local_scalar_t__ sigmaM;
            (void) sigmaM;  // dummy to suppress unused var warning
            if (jacobian__)
                sigmaM = in__.scalar_lb_constrain(0,lp__);
            else
                sigmaM = in__.scalar_lb_constrain(0);

            vector<vector<local_scalar_t__> > Xt;
            size_t dim_Xt_0__ = nSub;
            Xt.resize(dim_Xt_0__);
            for (size_t k_0__ = 0; k_0__ < dim_Xt_0__; ++k_0__) {
                size_t dim_Xt_1__ = nObs;
                Xt[k_0__].reserve(dim_Xt_1__);
                for (size_t k_1__ = 0; k_1__ < dim_Xt_1__; ++k_1__) {
                    if (jacobian__)
                        Xt[k_0__].push_back(in__.scalar_constrain(lp__));
                    else
                        Xt[k_0__].push_back(in__.scalar_constrain());
                }
            }


            // transformed parameters



            // validate transformed parameters

            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning

            // model body

            current_statement_begin__ = 88;
            lp_accum__.add(lognormal_log<propto__>(Cl, 0, sdDef));
            current_statement_begin__ = 89;
            lp_accum__.add(lognormal_log<propto__>(Ka, 0, sdDef));
            current_statement_begin__ = 90;
            lp_accum__.add(lognormal_log<propto__>(Ke, 0, sdDef));
            current_statement_begin__ = 91;
            lp_accum__.add(lognormal_log<propto__>(sigmaP, 0, sdDef));
            current_statement_begin__ = 92;
            lp_accum__.add(lognormal_log<propto__>(sigmaM, 0, sdDef));
            current_statement_begin__ = 96;
            for (int ii = 1; ii <= nSub; ++ii) {

                current_statement_begin__ = 97;
                lp_accum__.add(PK1_SDE_lpdf<propto__>(get_base1(Xt,ii,"Xt",1), Ka, Ke, Cl, sigmaP, get_base1(D,ii,"D",1), get_base1(t,ii,"t",1), get_base1(dt,ii,"dt",1), pstream__));
            }
            current_statement_begin__ = 100;
            lp_accum__.add(normal_log<propto__>(to_array_1d(Yt), to_array_1d(Xt), sigmaM));

        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }

        lp_accum__.add(lp__);
        return lp_accum__.sum();

    } // log_prob()

    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }


    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("Cl");
        names__.push_back("Ka");
        names__.push_back("Ke");
        names__.push_back("sigmaP");
        names__.push_back("sigmaM");
        names__.push_back("Xt");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(nSub);
        dims__.push_back(nObs);
        dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;

        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);
        static const char* function__ = "model_PK1_Fixed_SDE_Noise_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double Cl = in__.scalar_lb_constrain(0);
        double Ka = in__.scalar_lb_constrain(0);
        double Ke = in__.scalar_lb_constrain(0);
        double sigmaP = in__.scalar_lb_constrain(0);
        double sigmaM = in__.scalar_lb_constrain(0);
        vector<vector<double> > Xt;
        size_t dim_Xt_0__ = nSub;
        Xt.resize(dim_Xt_0__);
        for (size_t k_0__ = 0; k_0__ < dim_Xt_0__; ++k_0__) {
            size_t dim_Xt_1__ = nObs;
            for (size_t k_1__ = 0; k_1__ < dim_Xt_1__; ++k_1__) {
                Xt[k_0__].push_back(in__.scalar_constrain());
            }
        }
        vars__.push_back(Cl);
        vars__.push_back(Ka);
        vars__.push_back(Ke);
        vars__.push_back(sigmaP);
        vars__.push_back(sigmaM);
            for (int k_1__ = 0; k_1__ < nObs; ++k_1__) {
                for (int k_0__ = 0; k_0__ < nSub; ++k_0__) {
                vars__.push_back(Xt[k_0__][k_1__]);
                }
            }

        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {



            // validate transformed parameters

            // write transformed parameters
            if (include_tparams__) {
            }
            if (!include_gqs__) return;
            // declare and define generated quantities



            // validate generated quantities

            // write generated quantities
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    static std::string model_name() {
        return "model_PK1_Fixed_SDE_Noise";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "Cl";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "Ka";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "Ke";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigmaP";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigmaM";
        param_names__.push_back(param_name_stream__.str());
        for (int k_1__ = 1; k_1__ <= nObs; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= nSub; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "Xt" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
        }


        if (!include_gqs__) return;
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "Cl";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "Ka";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "Ke";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigmaP";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigmaM";
        param_names__.push_back(param_name_stream__.str());
        for (int k_1__ = 1; k_1__ <= nObs; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= nSub; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "Xt" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
        }


        if (!include_gqs__) return;
    }

}; // model

}

typedef model_PK1_Fixed_SDE_Noise_namespace::model_PK1_Fixed_SDE_Noise stan_model;


#endif