#include <TMB.hpp>
using namespace density;

//////////// helper fun: find NA locations in vector ///////////
template <class Type>
vector<Type> is_not_na(vector<Type> x){
	vector<Type> y(x.size());
	y.fill(Type(1.0));
    for(int i=0; i<x.size(); i++){
		if( R_IsNA(asDouble(x(i))) ){
			y(i) = Type(0.0);
		}
	}
	return y;
}

//////////// helper fun: extract non-NAs from vector ///////////
template<class Type>
vector<Type> remove_nas(vector<Type> data_vector, int number_of_datas, vector<Type> na_bool){
  int ii = 0;
	vector<Type> y_reduced(number_of_datas);
	for(int i=0; i < data_vector.size(); i++){
		if(na_bool(i) == Type(1.0)){
			y_reduced(ii) = data_vector(i);
			ii++;
		}
	}
  return y_reduced;
}

//////////// helper fun: construct permutation matrix ///////////
template <class Type>
matrix<Type> construct_permutation_matrix(int number_of_datas, int m, vector<Type> na_bool){
	matrix<Type> E(number_of_datas,m);
	E.setZero();
	/**/
	int j=0;
	for(int i=0; i < m; i++){
		if(na_bool(i) == Type(1.0)){ /*if p(i) is 1 then include by setting 1 in diagonal of matrix*/
			E(j,i) = Type(1.0);
			j += 1;
		}
	}
	return E;
}

//////////// Loss function ///////////
template<class Type>
Type lossfunction__(Type x, vector<Type> tukeypars, Type huber_c, int lossFunc){
  Type loss;
  if(lossFunc==1){
    Type a = tukeypars(0);
    Type b = tukeypars(1);
    Type c = tukeypars(2);
    Type d = tukeypars(3);
    loss = d * ( (Type(1.0)/(Type(1.0)+exp(-a*(x-b)))) + c );
  } else if (lossFunc==2){
    Type c_squared = pow(huber_c,2);
    loss = c_squared * (sqrt(1 + (x / c_squared)) - 1);
  } else {
    loss = x;
  }
  return(loss);
}

//////////// Map estimation helper ///////////
template<class Type>
vector<Type> get_free_pars__(vector<int> mapints, int sum_mapints, vector<Type> parvec) {
	vector<Type> ans(sum_mapints);
	int j=0;
	for(int i=0;i<mapints.size();i++){
		if(mapints(i)==1){
			ans(j) = parvec(i);
			j += 1;
		}
	}
	return(ans);
}
template<class Type>
vector<Type> f__(Type logtheta, Type mu, Type x){
	vector<Type> ans(1);
	ans(0) = exp(logtheta) * (mu - x);
	return ans;
}
template<class Type>
matrix<Type> dfdx__(Type logtheta){
	matrix<Type> ans(1,1);
	ans(0,0) = -exp(logtheta);
	return ans;
}
template<class Type>
matrix<Type> g__(Type logsigma_x){
	matrix<Type> ans(1,1);
	ans(0,0) = exp(logsigma_x);
	return ans;
}
template<class Type>
vector<Type> h__(Type x){
	vector<Type> ans(1);
	ans(0) = x;
	return ans;
}
template<class Type>
matrix<Type> dhdx__(Type a){
	matrix<Type> ans(1,1);
	ans(0,0) = 1;
	return ans;
}
template<class Type>
matrix<Type> obsvarFun__(Type logsigma_y){
	matrix<Type> V(1,1);
	V(0,0) = pow(exp(logsigma_y), 2);
	return V;
}

//////////// objective function ///////////
template<class Type>
Type objective_function<Type>::operator() ()
{

//// observations ////
	 DATA_VECTOR(y);

//// inputs ////
	 DATA_VECTOR(t);

//// initial state ////
	 DATA_VECTOR(X0__);
	 DATA_MATRIX(P0__);
	 DATA_VECTOR(dt__);
	 DATA_IVECTOR(N__);

//// loss parameters ////
	 DATA_VECTOR(tukey_pars__);
	 DATA_INTEGER(which_loss__);
	 DATA_SCALAR(loss_c_value__);

//// map estimation ////
	 DATA_INTEGER(map_bool__);

//// parameters ////
	 PARAMETER(logtheta);
	 PARAMETER(mu);
	 PARAMETER(logsigma_x);
	 PARAMETER(logsigma_y);

//// constants ////

//// system size ////
	 DATA_INTEGER(n__);
	 DATA_INTEGER(m__);

//////////// storage variables ///////////
	 vector<vector<Type>> xPrior(t.size());
	 vector<matrix<Type>> pPrior(t.size());
	 vector<vector<Type>> xPost(t.size());
	 vector<matrix<Type>> pPost(t.size());
	 vector<vector<Type>> Innovation(t.size());
	 vector<matrix<Type>> InnovationCovariance(t.size());

//////////// set initial value ///////////
	 vector<Type> x0__ = X0__;
	 matrix<Type> p0__ = P0__;
	 xPrior(0) = X0__;
	 xPost(0) = X0__;
	 pPrior(0) = P0__;
	 pPost(0) = P0__;

	 //////////// initialize variables ///////////
	 int s__;
	 Type half_log2PI = Type(0.5)*log(2*M_PI);
	 Type nll__ = 0;
	 vector<Type> data_vector__(m__),na_bool__,e__,y__,F__,H__;
	 matrix<Type> C__,R__,K__,E__,V__,Ri__,A__,G__;

	 //////////// identity matrix ///////////
	 matrix<Type> I__(n__,n__);
	 I__.setIdentity();

	 //////////// MAIN LOOP OVER TIME POINTS ///////////
	 for(int i=0 ; i<t.size()-1 ; i++){

		 //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
		 for(int j=0 ; j<N__(i) ; j++){
			 F__  = f__(logtheta, mu, x0__(0));
			 A__  = dfdx__(logtheta);
			 G__  = g__(logsigma_x);
			 x0__ = x0__ + F__ * dt__(i);
			 p0__ = p0__ + ( A__*p0__ + p0__*A__.transpose() + G__*G__.transpose() ) * dt__(i);
		 }
		 xPrior(i+1) = x0__;
		 pPrior(i+1) = p0__;

		 //////////// DATA-UPDATE ///////////
		 data_vector__ << y(i+1);
		 na_bool__ = is_not_na(data_vector__);
		 s__ = CppAD::Integer(sum(na_bool__));
		 if( s__ > 0 ){
			 y__  = remove_nas(data_vector__, s__, na_bool__);
			 E__  = construct_permutation_matrix(s__, m__, na_bool__);
			 H__  = h__(x0__(0));
			 C__  = E__ * dhdx__(Type(0.0));
			 e__  = y__ - E__ * H__;
			 V__  = E__ * obsvarFun__(logsigma_y) * E__.transpose();
			 R__  = C__ * p0__ * C__.transpose() + V__;
			 Ri__ = R__.inverse();
			 K__ 	= p0__ * C__.transpose() * Ri__;
			 x0__ = x0__ + K__*e__;
			 p0__ = (I__ - K__ * C__) * p0__ * (I__ - K__ * C__).transpose() + K__* V__ * K__.transpose();
			 nll__ += Type(0.5)*atomic::logdet(R__) + Type(0.5)*lossfunction__((e__*(Ri__*e__)).sum(),tukey_pars__,loss_c_value__,which_loss__) + half_log2PI * asDouble(s__);
			 Innovation(i+1) = e__;
			 InnovationCovariance(i+1) = R__;
		 }
		 xPost(i+1) = x0__;
		 pPost(i+1) = p0__;
	 }

		 //////////// MAP CONTRIBUTION ///////////
	 if(map_bool__==1){
		 DATA_VECTOR(map_mean__);
		 DATA_MATRIX(map_cov__);
		 DATA_IVECTOR(map_ints__);
		 DATA_INTEGER(sum_map_ints__);
		 vector<Type> parvec__(4);
		 vector<Type> map_pars__;
		 parvec__ << logtheta, mu, logsigma_x, logsigma_y;
		 map_pars__ = get_free_pars__(map_ints__,sum_map_ints__,parvec__);
		 vector<Type> pars_eps__ = map_pars__ - map_mean__;
		 matrix<Type> map_invcov__ = map_cov__.inverse();
		 Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();
		 nll__ += map_nll__;
		 REPORT(map_nll__);
		 REPORT(map_pars__);
		 REPORT(pars_eps__);
	 }

	 //////////// Return/Report //////////////
	 REPORT(Innovation);
	 REPORT(InnovationCovariance);
	 REPORT(xPrior);
	 REPORT(xPost);
	 REPORT(pPrior);
	 REPORT(pPost);
	 return nll__;
}
