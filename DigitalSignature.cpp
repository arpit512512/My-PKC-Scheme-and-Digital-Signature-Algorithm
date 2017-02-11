#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pari/pari.h>
//#include <pari.h>
#include <time.h>
#include <string.h>
#define PARI_OLD_NAMES
#include <iostream>

using namespace std;

// Some functions for random sample generation
double Uniform(void) {
	return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
} 
/*********************************
Standard Normal Distribution: use Box-Muller method
**********************************/
double Normal(void) {
	return sqrt( -log(Uniform())*2.0 ) * sin( 2.0*M_PI*Uniform() );
}
/*******************************
Gaussian Distribution
********************************/
double Gauss(double mu, double sigma) {   
	double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
	return mu + sigma*z;
}
/********************************************************************
Error sampling function
GEN Sample(
	int n			: polynomial dimension
	double sigma	: standard deviation parameter 

Contents :  Create n dimension random number vector(signed integer) based on Gaussian distribution (0, σ)
********************************************************************/
GEN Sample(int n, double sigma)
{
	GEN ret	= cgetg(n + 1, t_VEC);
	double z;
	int i;
	
	for (i = 1; i <= n; i++) {
		z = Gauss(0, sigma); z = abs(round(z)); /*absolute value of Gaussian distribution */
		ret[i] = (long) stoi((long) z);
	}
	
	return ret;
}

GEN randomElement(int n){
	GEN ret;
	ret = cgetg(n + 1, t_VEC);
	for(int i=0; i<n; i++){
		gel(ret, i+1) = stoi(rand());
	}
	//pari_printf("ret: %s\n", GENtostr(ret));
	return ret;
}

GEN randomElementredbyt(int n, GEN t){
	GEN ret;
	ret = cgetg(n + 1, t_VEC);
	for(int i=0; i<n; i++){
		gel(ret, i+1) = lift(gmodulo(stoi(rand()), t));
	}
	//pari_printf("ret: %s\n", GENtostr(ret));
	return ret;
}


int main(){
	pari_init(3000000000,2);
	//GEN a;
	//a = gadd(stoi(45), stoi(1));
	//pari_printf("%s\n", GENtostr(a));
	//cout<<rand()<<" "<<rand()<<endl;
	GEN n, q, t, bit, a, b, c, d, e, f, p;
	bit = powii(stoi(2), stoi(215));
	cout<<GENtostr(bit)<<endl;
	a = nextprime(bit);
	b = nextprime(gadd(a, stoi(1)));
	GEN temp, temp1, temp2;
	temp = gmul(a, b);
	temp = gmul(stoi(2), temp);
	//cout<<isprime(stoi(7))<<endl;
	while(!isprime(gadd(stoi(1), temp))){
		//cout<<GENtostr(a)<<endl;
		a = nextprime(gadd(b, stoi(1)));
		b = nextprime(gadd(a, stoi(1)));
		temp = gmul(a, b);
		temp = gmul(stoi(2), temp);
	}
	e = gadd(stoi(1), temp);
	GEN eminus = gsub(e, stoi(1));
	cout<<GENtostr(e)<<endl;
	temp = stoi(2001);
	p = e;
	cout<<endl<<GENtostr(Fp_pow(temp, eminus, p));
	cout<<endl<<GENtostr(Fp_pow(temp, a, p));
	cout<<endl<<GENtostr(Fp_pow(temp, b, p));
	cout<<endl<<GENtostr(Fp_pow(temp, gmul(b, a), p));
	GEN a1, alpha, alphaa, m, k, gamma, delta, beta, l, pminus;
	pminus = gsub(p, stoi(1));
	alpha = temp;
	cout<<"a\n";
	l = stoi(rand()+1);
	while(1){
		if(itos(lift(gmodulo(l, stoi(2))))==1)
			break;
		l = stoi(rand()+1);
	}
	k = stoi(rand()+1);
	while(1){
		if(gcmp(gcdii(k, pminus), stoi(1))==0)
			break;
		k = stoi(rand()+1);
	}
	cout<<"a\n";
	
	beta = Fp_pow(alpha, l, p);
	GEN r, kinv, s, exp, alphainv, temp11;
	cout<<"a\n";
	r = Fp_pow(alpha, k, p);
	kinv = ginvmod(k, pminus);
	cout<<"b\n";
	temp11 = gmul(l, r);
	//temp1 = lift(gmodulo(temp1, pminus));
	
	cout<<"a\n";
	clock_t start, end;
	int noofmessages = 200000;
	start = clock();
	
	for(int i=0; i<noofmessages; i++){
		m = stoi(rand()+1);
		//cout<<itos(m)<<endl;
		c = stoi(rand()+1);
		//cout<<itos(c)<<endl;
	while(1){
	if(gcmp(k, c)!=0 && itos(c)%2==1)
		break;
	c = stoi(rand()+1);
	}
	while(1){
		if(itos(m)%2==1)
		break;
		m = stoi(rand()+1);
	}
	//cout<<"b\n";
		temp = gmul(c, m);
		temp2 = temp;
		//temp = lift(gmodulo(temp, pminus));
		temp = gsub(temp, temp11);
		temp = lift(gmodulo(temp, pminus));
		temp1 = temp;
		//cout<<GENtostr(temp1)<<endl;
		//
		s = gmul(temp, kinv);
		//cout<<gcmp(lift(gmodulo(gmul(s, k), pminus)), lift(gmodulo(temp1, pminus)))<<endl;
	//s = lift(gmodulo(s, pminus));
	temp = temp2;
	alphainv = ginvmod(alpha, p);
	//cout<<"b\n";
	//cout<<GENtostr(gmodulo(gmul(Fp_pow(alphainv, temp, p), Fp_pow(alpha, temp, p)), p))<<endl;
	//cout<<"a-1 -> "<<GENtostr(alphainv)<<endl;
	exp = Fp_pow(alphainv, temp, p);
	//cout<<gcmp(lift(gmodulo(gmul(Fp_pow(r, s, p), Fp_pow(beta, r, p)), p)), Fp_pow(alpha, temp2, p))<<endl;
	//cout<<GENtostr(gcdii(gmul(c, m), pminus))<<endl;
	
	//cout<<GENtostr(lift(gmodulo(gmul(exp, lift(gmodulo(gmul(alpha, temp), p))), p)))<<endl;
	temp = Fp_pow(r, s, p);
	temp1 = Fp_pow(beta, r, p);
	temp = gmul(temp, temp1);
	temp = lift(gmodulo(temp, p));
	temp = gmul(temp, exp);
	temp = lift(gmodulo(temp, p));
	//cout<<"Answer -> "<<GENtostr(temp)<<endl;
	}
	
	//a1 = a;
	//alphaa = Fp_pow(alpha, a1, p);
	/*beta = Fp_pow(alpha, l, p);
	clock_t start, end;
	int noofmessages = 20000;
	start = clock();
	k = powii(stoi(2), stoi(250));
	cout<<bit_prec(k)<<endl;
	cout<<GENtostr(k)<<endl;
	for(int i=0; i<noofmessages; i++){
		m = stoi(rand());
	
	gamma = Fp_pow(alpha, k, p);
	delta = Fp_pow(alphaa, k, p);
	delta = gmul(delta, m);
	delta = lift(gmodulo(delta, p));
	temp = ginvmod(gamma, p);
	temp1 = Fp_pow(temp, a1, p);
	//temp = gsub(p, stoi(1));
	//temp = gsub(temp, a);
	//temp2 = Fp_pow(gamma, temp, p);
	temp2 = gmul(temp1, delta);
	temp2 = lift(gmodulo(temp2, p));
	cout<<GENtostr(temp2)<<endl;
	}
	
	/*c = nextprime(gadd(b, stoi(1)));
	d = nextprime(gadd(c, stoi(1)));
	temp = gmul(c, d);
	temp = gmul(stoi(2), temp);
	while(!isprime(gadd(stoi(1), temp))){
		//cout<<GENtostr(a)<<endl;
		c = nextprime(gadd(d, stoi(1)));
		d = nextprime(gadd(c, stoi(1)));
		temp = gmul(c, d);
		temp = gmul(stoi(2), temp);
	}
	f = gadd(stoi(1), temp);
	cout<<GENtostr(f)<<"\n";
	cout<<isprime(e)<<"    "<<isprime(f)<<endl;
	GEN eminus, elorder1, elorder2, fminus;
	eminus = gsub(e, stoi(1));
	fminus = gsub(f, stoi(1));
	n = gmul(e, f);
	GEN g = stoi(200000);
	temp = gsub(f, stoi(1));
	temp1 = Fp_pow(g, temp, n);
	cout<<GENtostr(temp1)<<endl;
	elorder1 = temp1;
	cout<<GENtostr(Fp_pow(temp1, eminus, n))<<endl;
	g = gadd(g, stoi(1));
	temp = gmul(stoi(2), fminus);
	temp = gmul(temp, a);
	while(1){
		temp1 = Fp_pow(g, temp, n);
		elorder2 = temp1;
		cout<<GENtostr(temp1)<<endl;
		temp2 = gmul(a,b);
		temp2 = gmul(temp2, stoi(2));
		if(gcmp(Fp_pow(temp1, b, n), stoi(1))==0)
			break;
		g = gadd(g, stoi(1));
	}
	cout<<GENtostr(g)<<endl;
	GEN exp, c1, alpha, beta, m, k;
	
	k = stoi(5);
	alpha = elorder1;
	beta = Fp_pow(alpha, gadd(b, stoi(1)), n);
	GEN betak, alpha2k;
	betak = Fp_pow(beta, k, n);
	alpha2k = Fp_pow(alpha, gmul(stoi(2), k), n);
	clock_t start, end;
	cout<<GENtostr(beta)<<endl;
	cout<<"Looping\n";
	start = clock();
	int noofmessages = 10000;
	for(int i=0; i<noofmessages; i++){
		c1 = stoi(rand());
	m = stoi(rand());
	exp = Fp_pow(elorder2, c1, n);
	GEN delta, gamma;
	delta = betak;
	delta = gmul(delta, exp);
	delta = gmul(delta, m);
	gamma = gmul(exp, exp);
	temp = alpha2k;
	gamma = gmul(gamma, temp);
	GEN gammainv;
	//cout<<"Finding inverse\n";
	gammainv = ginvmod(gamma, n);
	//cout<<GENtostr(gmodulo(gmul(gammainv, gamma), n))<<endl;
	temp = gmul(delta, delta);
	temp = gmul(temp, gammainv);
	temp = lift(gmodulo(temp, n));
	//cout<<GENtostr(n)<<endl<<"Answer\n";
	temp1 = Fp_sqrt(temp, e);
	temp2 = Fp_sqrt(temp, f);
	GEN fac1, fac2, fac;
	fac1 = cgetg(3, t_COL);
	fac2 = cgetg(3, t_COL);
	gel(fac1, 1) = e;
	gel(fac1, 2) = f;
	gel(fac2, 1) = stoi(1);
	gel(fac2, 2) = stoi(1);
	fac = cgetg(3, t_MAT);
	gel(fac, 1) = fac1;
	gel(fac, 2) = fac2;
	temp = Zn_sqrt(temp, fac);
	//cout<<GENtostr(temp)<<endl<<GENtostr(temp1)<<endl<<GENtostr(temp2)<<endl;
	//cout<<GENtostr(gmodulo(e, stoi(4)))<<"    "<<GENtostr(gmodulo(f, stoi(4)))<<endl;
	}*/
	end = clock();
	int timetaken = ((end - start)/CLOCKS_PER_SEC);
	cout<<endl<<timetaken<<" seconds\n"<<(float)timetaken/(float)noofmessages<<" seconds per message"<<endl;
	//cout<<bit_prec(p)<<endl;
	//while(1){
	//	if(gcmp(Fp_pow(g, gmul(eminus, gmul(c,d)), n), stoi(1))==0)
	//		break;
	//	g = gadd(g, stoi(1));
	//	cout<<GENtostr(g)<<endl;
	//}
	//cout<<GENtostr(g)<<endl;
	pari_close();
	return 0;
}