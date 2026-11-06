#define main softframe_demo_main
#include "../softframe/main.cpp"
#undef main
int main(){
 for(unsigned k=0;k<1024;k++){
  RNG r{123};r.setup_quasi();r.sample_index(k);
  for(unsigned d=0;d<32;d++){
   unsigned bits=k^(k>>1),value=r.shifts[d];
   for(unsigned b=0;b<32;b++)if(bits&(1u<<b))value^=sobol_directions[d][b];
   assert(r.get()==(double(value)+.5)*0x1.0p-32);
  }
 }
 Scene scene;RNG rng{7654321};
 for(int i=0;i<500;i++){V a{rng.get()*10-5,rng.get()*10-5,rng.get()*10};scene.tris.push_back({a,a+V{.4,0,0},a+V{0,.4,.02}});}
 scene.build(0,scene.tris.size());scene.accelerate();
 for(int i=0;i<10000;i++){Ray ray{{0,0,-2},norm({rng.get()*2-1,rng.get()*2-1,1})};double x=1e30,y=1e30;V a,b;bool fast=scene.trace(ray,x,a),ref=scene.trace(ray,y,b,true);assert(fast==ref);if(fast)assert(std::abs(x-y)<1e-4);}
 Scene single;single.tris.push_back({{-1,-1,1},{1,-1,1},{0,1,1}});single.build(0,1);single.accelerate();
 for(Ray ray:std::vector<Ray>{{{0,0,0},{0,0,1}},{{0,0,0},{0,0,-1}},{{0,0,0},{1,0,0}},{{-1,-1,0},{0,0,1}},{{0,0,1},{0,0,1}}}){double x=1e30,y=1e30;V a,b;assert(single.trace(ray,x,a)==single.trace(ray,y,b,true));assert(std::abs(x-y)<1e-4);}
 std::cout<<"{\"random_reference_checks\":10000,\"axis_boundary_checks\":5,\"sobol_exact_value_checks\":32768,\"passed\":true}\n";
}
