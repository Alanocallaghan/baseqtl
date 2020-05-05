// Generated by rstantools.  Do not edit by hand.

/*
    baseqtl is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    baseqtl is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with baseqtl.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.18.0
#include <stan/model/model_header.hpp>
namespace model_GT_nb_ase_refbias_namespace {
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
    reader.add_event(0, 0, "start", "model_GT_nb_ase_refbias");
    reader.add_event(96, 94, "end", "model_GT_nb_ase_refbias");
    return reader;
}
#include <stan_meta_header.hpp>
class model_GT_nb_ase_refbias : public prob_grad {
private:
    int N;
    int A;
    int L;
    int K;
    int k;
    vector<int> Y;
    vector<int> g;
    vector<int> gase;
    vector<int> m;
    vector<int> n;
    vector_d pH;
    vector_d ai0;
    vector_d sdai0;
    vector<int> s;
    matrix_d cov;
    vector_d aveP;
    vector_d sdP;
    vector_d mixP;
public:
    model_GT_nb_ase_refbias(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_GT_nb_ase_refbias(stan::io::var_context& context__,
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
        static const char* function__ = "model_GT_nb_ase_refbias_namespace::model_GT_nb_ase_refbias";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        // initialize member variables
        try {
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "A", "int", context__.to_vec());
            A = int(0);
            vals_i__ = context__.vals_i("A");
            pos__ = 0;
            A = vals_i__[pos__++];
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "L", "int", context__.to_vec());
            L = int(0);
            vals_i__ = context__.vals_i("L");
            pos__ = 0;
            L = vals_i__[pos__++];
            current_statement_begin__ = 7;
            context__.validate_dims("data initialization", "K", "int", context__.to_vec());
            K = int(0);
            vals_i__ = context__.vals_i("K");
            pos__ = 0;
            K = vals_i__[pos__++];
            current_statement_begin__ = 8;
            context__.validate_dims("data initialization", "k", "int", context__.to_vec());
            k = int(0);
            vals_i__ = context__.vals_i("k");
            pos__ = 0;
            k = vals_i__[pos__++];
            current_statement_begin__ = 9;
            validate_non_negative_index("Y", "N", N);
            context__.validate_dims("data initialization", "Y", "int", context__.to_vec(N));
            validate_non_negative_index("Y", "N", N);
            Y = std::vector<int>(N,int(0));
            vals_i__ = context__.vals_i("Y");
            pos__ = 0;
            size_t Y_limit_0__ = N;
            for (size_t i_0__ = 0; i_0__ < Y_limit_0__; ++i_0__) {
                Y[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 10;
            validate_non_negative_index("g", "N", N);
            context__.validate_dims("data initialization", "g", "int", context__.to_vec(N));
            validate_non_negative_index("g", "N", N);
            g = std::vector<int>(N,int(0));
            vals_i__ = context__.vals_i("g");
            pos__ = 0;
            size_t g_limit_0__ = N;
            for (size_t i_0__ = 0; i_0__ < g_limit_0__; ++i_0__) {
                g[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 11;
            validate_non_negative_index("gase", "A", A);
            context__.validate_dims("data initialization", "gase", "int", context__.to_vec(A));
            validate_non_negative_index("gase", "A", A);
            gase = std::vector<int>(A,int(0));
            vals_i__ = context__.vals_i("gase");
            pos__ = 0;
            size_t gase_limit_0__ = A;
            for (size_t i_0__ = 0; i_0__ < gase_limit_0__; ++i_0__) {
                gase[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 12;
            validate_non_negative_index("m", "A", A);
            context__.validate_dims("data initialization", "m", "int", context__.to_vec(A));
            validate_non_negative_index("m", "A", A);
            m = std::vector<int>(A,int(0));
            vals_i__ = context__.vals_i("m");
            pos__ = 0;
            size_t m_limit_0__ = A;
            for (size_t i_0__ = 0; i_0__ < m_limit_0__; ++i_0__) {
                m[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 13;
            validate_non_negative_index("n", "L", L);
            context__.validate_dims("data initialization", "n", "int", context__.to_vec(L));
            validate_non_negative_index("n", "L", L);
            n = std::vector<int>(L,int(0));
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            size_t n_limit_0__ = L;
            for (size_t i_0__ = 0; i_0__ < n_limit_0__; ++i_0__) {
                n[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 14;
            validate_non_negative_index("pH", "L", L);
            context__.validate_dims("data initialization", "pH", "vector_d", context__.to_vec(L));
            validate_non_negative_index("pH", "L", L);
            pH = vector_d(static_cast<Eigen::VectorXd::Index>(L));
            vals_r__ = context__.vals_r("pH");
            pos__ = 0;
            size_t pH_i_vec_lim__ = L;
            for (size_t i_vec__ = 0; i_vec__ < pH_i_vec_lim__; ++i_vec__) {
                pH[i_vec__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 15;
            validate_non_negative_index("ai0", "L", L);
            context__.validate_dims("data initialization", "ai0", "vector_d", context__.to_vec(L));
            validate_non_negative_index("ai0", "L", L);
            ai0 = vector_d(static_cast<Eigen::VectorXd::Index>(L));
            vals_r__ = context__.vals_r("ai0");
            pos__ = 0;
            size_t ai0_i_vec_lim__ = L;
            for (size_t i_vec__ = 0; i_vec__ < ai0_i_vec_lim__; ++i_vec__) {
                ai0[i_vec__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 16;
            validate_non_negative_index("sdai0", "L", L);
            context__.validate_dims("data initialization", "sdai0", "vector_d", context__.to_vec(L));
            validate_non_negative_index("sdai0", "L", L);
            sdai0 = vector_d(static_cast<Eigen::VectorXd::Index>(L));
            vals_r__ = context__.vals_r("sdai0");
            pos__ = 0;
            size_t sdai0_i_vec_lim__ = L;
            for (size_t i_vec__ = 0; i_vec__ < sdai0_i_vec_lim__; ++i_vec__) {
                sdai0[i_vec__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 17;
            validate_non_negative_index("s", "A", A);
            context__.validate_dims("data initialization", "s", "int", context__.to_vec(A));
            validate_non_negative_index("s", "A", A);
            s = std::vector<int>(A,int(0));
            vals_i__ = context__.vals_i("s");
            pos__ = 0;
            size_t s_limit_0__ = A;
            for (size_t i_0__ = 0; i_0__ < s_limit_0__; ++i_0__) {
                s[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 18;
            validate_non_negative_index("cov", "N", N);
            validate_non_negative_index("cov", "(1 + K)", (1 + K));
            context__.validate_dims("data initialization", "cov", "matrix_d", context__.to_vec(N,(1 + K)));
            validate_non_negative_index("cov", "N", N);
            validate_non_negative_index("cov", "(1 + K)", (1 + K));
            cov = matrix_d(static_cast<Eigen::VectorXd::Index>(N),static_cast<Eigen::VectorXd::Index>((1 + K)));
            vals_r__ = context__.vals_r("cov");
            pos__ = 0;
            size_t cov_m_mat_lim__ = N;
            size_t cov_n_mat_lim__ = (1 + K);
            for (size_t n_mat__ = 0; n_mat__ < cov_n_mat_lim__; ++n_mat__) {
                for (size_t m_mat__ = 0; m_mat__ < cov_m_mat_lim__; ++m_mat__) {
                    cov(m_mat__,n_mat__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 19;
            validate_non_negative_index("aveP", "k", k);
            context__.validate_dims("data initialization", "aveP", "vector_d", context__.to_vec(k));
            validate_non_negative_index("aveP", "k", k);
            aveP = vector_d(static_cast<Eigen::VectorXd::Index>(k));
            vals_r__ = context__.vals_r("aveP");
            pos__ = 0;
            size_t aveP_i_vec_lim__ = k;
            for (size_t i_vec__ = 0; i_vec__ < aveP_i_vec_lim__; ++i_vec__) {
                aveP[i_vec__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 20;
            validate_non_negative_index("sdP", "k", k);
            context__.validate_dims("data initialization", "sdP", "vector_d", context__.to_vec(k));
            validate_non_negative_index("sdP", "k", k);
            sdP = vector_d(static_cast<Eigen::VectorXd::Index>(k));
            vals_r__ = context__.vals_r("sdP");
            pos__ = 0;
            size_t sdP_i_vec_lim__ = k;
            for (size_t i_vec__ = 0; i_vec__ < sdP_i_vec_lim__; ++i_vec__) {
                sdP[i_vec__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 21;
            validate_non_negative_index("mixP", "k", k);
            context__.validate_dims("data initialization", "mixP", "vector_d", context__.to_vec(k));
            validate_non_negative_index("mixP", "k", k);
            mixP = vector_d(static_cast<Eigen::VectorXd::Index>(k));
            vals_r__ = context__.vals_r("mixP");
            pos__ = 0;
            size_t mixP_i_vec_lim__ = k;
            for (size_t i_vec__ = 0; i_vec__ < mixP_i_vec_lim__; ++i_vec__) {
                mixP[i_vec__] = vals_r__[pos__++];
            }
            // validate, data variables
            current_statement_begin__ = 4;
            check_greater_or_equal(function__,"N",N,0);
            current_statement_begin__ = 5;
            check_greater_or_equal(function__,"A",A,0);
            current_statement_begin__ = 6;
            check_greater_or_equal(function__,"L",L,0);
            current_statement_begin__ = 7;
            check_greater_or_equal(function__,"K",K,0);
            current_statement_begin__ = 8;
            check_greater_or_equal(function__,"k",k,0);
            current_statement_begin__ = 9;
            current_statement_begin__ = 10;
            current_statement_begin__ = 11;
            current_statement_begin__ = 12;
            current_statement_begin__ = 13;
            current_statement_begin__ = 14;
            current_statement_begin__ = 15;
            current_statement_begin__ = 16;
            current_statement_begin__ = 17;
            current_statement_begin__ = 18;
            current_statement_begin__ = 19;
            current_statement_begin__ = 20;
            current_statement_begin__ = 21;
            // initialize data variables
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 28;
            validate_non_negative_index("betas", "K", K);
            num_params_r__ += K;
            current_statement_begin__ = 29;
            ++num_params_r__;
            current_statement_begin__ = 30;
            ++num_params_r__;
            current_statement_begin__ = 31;
            ++num_params_r__;
            current_statement_begin__ = 32;
            validate_non_negative_index("rai0", "L", L);
            num_params_r__ += L;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_GT_nb_ase_refbias() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        if (!(context__.contains_r("betas")))
            throw std::runtime_error("variable betas missing");
        vals_r__ = context__.vals_r("betas");
        pos__ = 0U;
        validate_non_negative_index("betas", "K", K);
        context__.validate_dims("initialization", "betas", "vector_d", context__.to_vec(K));
        vector_d betas(static_cast<Eigen::VectorXd::Index>(K));
        for (int j1__ = 0U; j1__ < K; ++j1__)
            betas(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_unconstrain(betas);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable betas: ") + e.what());
        }
        if (!(context__.contains_r("bj")))
            throw std::runtime_error("variable bj missing");
        vals_r__ = context__.vals_r("bj");
        pos__ = 0U;
        context__.validate_dims("initialization", "bj", "double", context__.to_vec());
        double bj(0);
        bj = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(bj);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable bj: ") + e.what());
        }
        if (!(context__.contains_r("phi")))
            throw std::runtime_error("variable phi missing");
        vals_r__ = context__.vals_r("phi");
        pos__ = 0U;
        context__.validate_dims("initialization", "phi", "double", context__.to_vec());
        double phi(0);
        phi = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,phi);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable phi: ") + e.what());
        }
        if (!(context__.contains_r("theta")))
            throw std::runtime_error("variable theta missing");
        vals_r__ = context__.vals_r("theta");
        pos__ = 0U;
        context__.validate_dims("initialization", "theta", "double", context__.to_vec());
        double theta(0);
        theta = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,theta);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable theta: ") + e.what());
        }
        if (!(context__.contains_r("rai0")))
            throw std::runtime_error("variable rai0 missing");
        vals_r__ = context__.vals_r("rai0");
        pos__ = 0U;
        validate_non_negative_index("rai0", "L", L);
        context__.validate_dims("initialization", "rai0", "vector_d", context__.to_vec(L));
        vector_d rai0(static_cast<Eigen::VectorXd::Index>(L));
        for (int j1__ = 0U; j1__ < L; ++j1__)
            rai0(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_unconstrain(rai0);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable rai0: ") + e.what());
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
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  betas;
            (void) betas;  // dummy to suppress unused var warning
            if (jacobian__)
                betas = in__.vector_constrain(K,lp__);
            else
                betas = in__.vector_constrain(K);
            local_scalar_t__ bj;
            (void) bj;  // dummy to suppress unused var warning
            if (jacobian__)
                bj = in__.scalar_constrain(lp__);
            else
                bj = in__.scalar_constrain();
            local_scalar_t__ phi;
            (void) phi;  // dummy to suppress unused var warning
            if (jacobian__)
                phi = in__.scalar_lb_constrain(0,lp__);
            else
                phi = in__.scalar_lb_constrain(0);
            local_scalar_t__ theta;
            (void) theta;  // dummy to suppress unused var warning
            if (jacobian__)
                theta = in__.scalar_lb_constrain(0,lp__);
            else
                theta = in__.scalar_lb_constrain(0);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  rai0;
            (void) rai0;  // dummy to suppress unused var warning
            if (jacobian__)
                rai0 = in__.vector_constrain(L,lp__);
            else
                rai0 = in__.vector_constrain(L);
            // transformed parameters
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // model body
            {
            current_statement_begin__ = 37;
            validate_non_negative_index("lmu", "N", N);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  lmu(static_cast<Eigen::VectorXd::Index>(N));
            (void) lmu;  // dummy to suppress unused var warning
            stan::math::initialize(lmu, DUMMY_VAR__);
            stan::math::fill(lmu,DUMMY_VAR__);
            current_statement_begin__ = 38;
            validate_non_negative_index("p", "L", L);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  p(static_cast<Eigen::VectorXd::Index>(L));
            (void) p;  // dummy to suppress unused var warning
            stan::math::initialize(p, DUMMY_VAR__);
            stan::math::fill(p,DUMMY_VAR__);
            current_statement_begin__ = 39;
            local_scalar_t__ ebj;
            (void) ebj;  // dummy to suppress unused var warning
            stan::math::initialize(ebj, DUMMY_VAR__);
            stan::math::fill(ebj,DUMMY_VAR__);
            current_statement_begin__ = 40;
            int pos(0);
            (void) pos;  // dummy to suppress unused var warning
            stan::math::fill(pos, std::numeric_limits<int>::min());
            current_statement_begin__ = 41;
            validate_non_negative_index("ltmp", "L", L);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  ltmp(static_cast<Eigen::VectorXd::Index>(L));
            (void) ltmp;  // dummy to suppress unused var warning
            stan::math::initialize(ltmp, DUMMY_VAR__);
            stan::math::fill(ltmp,DUMMY_VAR__);
            current_statement_begin__ = 42;
            validate_non_negative_index("esum", "L", L);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  esum(static_cast<Eigen::VectorXd::Index>(L));
            (void) esum;  // dummy to suppress unused var warning
            stan::math::initialize(esum, DUMMY_VAR__);
            stan::math::fill(esum,DUMMY_VAR__);
            current_statement_begin__ = 43;
            validate_non_negative_index("esum0", "L", L);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  esum0(static_cast<Eigen::VectorXd::Index>(L));
            (void) esum0;  // dummy to suppress unused var warning
            stan::math::initialize(esum0, DUMMY_VAR__);
            stan::math::fill(esum0,DUMMY_VAR__);
            current_statement_begin__ = 44;
            validate_non_negative_index("lps", "k", k);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  lps(static_cast<Eigen::VectorXd::Index>(k));
            (void) lps;  // dummy to suppress unused var warning
            stan::math::initialize(lps, DUMMY_VAR__);
            stan::math::fill(lps,DUMMY_VAR__);
            current_statement_begin__ = 47;
            lp_accum__.add(gamma_log<propto__>(theta, 1, 0.10000000000000001));
            current_statement_begin__ = 48;
            lp_accum__.add(gamma_log<propto__>(phi, 1, 0.10000000000000001));
            current_statement_begin__ = 51;
            for (int i = 1; i <= k; ++i) {
                current_statement_begin__ = 52;
                stan::model::assign(lps, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            (normal_log(bj,get_base1(aveP,i,"aveP",1),get_base1(sdP,i,"sdP",1)) + get_base1(mixP,i,"mixP",1)), 
                            "assigning variable lps");
            }
            current_statement_begin__ = 54;
            lp_accum__.add(log_sum_exp(lps));
            current_statement_begin__ = 57;
            lp_accum__.add(normal_log<propto__>(get_base1(betas,1,"betas",1), 6, 4));
            current_statement_begin__ = 58;
            for (int i = 2; i <= K; ++i) {
                current_statement_begin__ = 59;
                lp_accum__.add(cauchy_log<propto__>(get_base1(betas,i,"betas",1), 0, 2.5));
            }
            current_statement_begin__ = 61;
            for (int i = 1; i <= L; ++i) {
                current_statement_begin__ = 62;
                lp_accum__.add(normal_log<propto__>(get_base1(rai0,i,"rai0",1), get_base1(ai0,i,"ai0",1), get_base1(sdai0,i,"sdai0",1)));
            }
            current_statement_begin__ = 66;
            stan::math::assign(ebj, stan::math::exp(bj));
            current_statement_begin__ = 68;
            stan::math::assign(lmu, multiply(stan::model::rvalue(cov, stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_min_max(2, cols(cov)), stan::model::nil_index_list())), "cov"),betas));
            current_statement_begin__ = 69;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 70;
                stan::model::assign(lmu, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            stan::model::deep_copy((logical_eq(stan::math::fabs(get_base1(g,i,"g",1)),1) ? stan::math::promote_scalar<local_scalar_t__>(((get_base1(lmu,i,"lmu",1) + stan::math::log1p(ebj)) - stan::math::log(2))) : stan::math::promote_scalar<local_scalar_t__>(get_base1(lmu,i,"lmu",1)) )), 
                            "assigning variable lmu");
                current_statement_begin__ = 71;
                stan::model::assign(lmu, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            stan::model::deep_copy((logical_eq(get_base1(g,i,"g",1),2) ? stan::math::promote_scalar<local_scalar_t__>((get_base1(lmu,i,"lmu",1) + bj)) : stan::math::promote_scalar<local_scalar_t__>(get_base1(lmu,i,"lmu",1)) )), 
                            "assigning variable lmu");
                current_statement_begin__ = 72;
                lp_accum__.add(neg_binomial_2_log(get_base1(Y,i,"Y",1),stan::math::exp(get_base1(lmu,i,"lmu",1)),phi));
            }
            current_statement_begin__ = 75;
            stan::math::assign(pos, 1);
            current_statement_begin__ = 76;
            stan::math::assign(esum, inv_logit(add(rai0,bj)));
            current_statement_begin__ = 77;
            stan::math::assign(esum0, inv_logit(rai0));
            current_statement_begin__ = 79;
            for (int i = 1; i <= A; ++i) {
                current_statement_begin__ = 80;
                for (int r = pos; r <= ((pos + get_base1(s,i,"s",1)) - 1); ++r) {
                    current_statement_begin__ = 82;
                    stan::model::assign(p, 
                                stan::model::cons_list(stan::model::index_uni(r), stan::model::nil_index_list()), 
                                (logical_eq(get_base1(gase,i,"gase",1),1) ? stan::math::promote_scalar<local_scalar_t__>(get_base1(esum,r,"esum",1)) : stan::math::promote_scalar<local_scalar_t__>(get_base1(esum0,r,"esum0",1)) ), 
                                "assigning variable p");
                    current_statement_begin__ = 83;
                    stan::model::assign(p, 
                                stan::model::cons_list(stan::model::index_uni(r), stan::model::nil_index_list()), 
                                stan::model::deep_copy((logical_eq(get_base1(gase,i,"gase",1),-(1)) ? stan::math::promote_scalar<local_scalar_t__>((1 - get_base1(esum,r,"esum",1))) : stan::math::promote_scalar<local_scalar_t__>(get_base1(p,r,"p",1)) )), 
                                "assigning variable p");
                    current_statement_begin__ = 85;
                    stan::model::assign(ltmp, 
                                stan::model::cons_list(stan::model::index_uni(r), stan::model::nil_index_list()), 
                                (beta_binomial_log(get_base1(n,r,"n",1),get_base1(m,i,"m",1),(get_base1(p,r,"p",1) * theta),((1 - get_base1(p,r,"p",1)) * theta)) + stan::math::log(get_base1(pH,r,"pH",1))), 
                                "assigning variable ltmp");
                }
                current_statement_begin__ = 87;
                lp_accum__.add(log_sum_exp(stan::model::rvalue(ltmp, stan::model::cons_list(stan::model::index_min_max(pos, ((pos + get_base1(s,i,"s",1)) - 1)), stan::model::nil_index_list()), "ltmp")));
                current_statement_begin__ = 88;
                stan::math::assign(pos, stan::model::deep_copy((pos + get_base1(s,i,"s",1))));
            }
            }
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
        names__.push_back("betas");
        names__.push_back("bj");
        names__.push_back("phi");
        names__.push_back("theta");
        names__.push_back("rai0");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(L);
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
        static const char* function__ = "model_GT_nb_ase_refbias_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        vector_d betas = in__.vector_constrain(K);
        double bj = in__.scalar_constrain();
        double phi = in__.scalar_lb_constrain(0);
        double theta = in__.scalar_lb_constrain(0);
        vector_d rai0 = in__.vector_constrain(L);
            for (int k_0__ = 0; k_0__ < K; ++k_0__) {
            vars__.push_back(betas[k_0__]);
            }
        vars__.push_back(bj);
        vars__.push_back(phi);
        vars__.push_back(theta);
            for (int k_0__ = 0; k_0__ < L; ++k_0__) {
            vars__.push_back(rai0[k_0__]);
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
        return "model_GT_nb_ase_refbias";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "betas" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "bj";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "phi";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "theta";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= L; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "rai0" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
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
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "betas" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "bj";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "phi";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "theta";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= L; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "rai0" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}
typedef model_GT_nb_ase_refbias_namespace::model_GT_nb_ase_refbias stan_model;
#endif
