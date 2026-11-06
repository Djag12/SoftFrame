#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>
#include <omp.h>
#ifdef SOFTFRAME_EMBREE
#include <embree4/rtcore.h>
#endif
constexpr double PI=3.141592653589793;
struct V {double x=0,y=0,z=0;V operator+(V b)const{return{x+b.x,y+b.y,z+b.z};}V operator-(V b)const{return{x-b.x,y-b.y,z-b.z};}V operator*(double s)const{return{x*s,y*s,z*s};}V operator*(V b)const{return{x*b.x,y*b.y,z*b.z};}double operator[](int i)const{return i==0?x:i==1?y:z;}};
double dot(V a,V b){return a.x*b.x+a.y*b.y+a.z*b.z;}V cross(V a,V b){return{a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}V norm(V a){return a*(1/std::sqrt(dot(a,a)));}
struct Ray{V o,d;};struct Tri{V a,b,c;V center()const{return(a+b+c)*(1./3);}};
struct Box{V lo{1e30,1e30,1e30},hi{-1e30,-1e30,-1e30};void add(V p){lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}bool hit(Ray r,double far)const{double near=1e-6;for(int k=0;k<3;k++){if(std::abs(r.d[k])<1e-14){if(r.o[k]<lo[k]||r.o[k]>hi[k])return false;continue;}double a=(lo[k]-r.o[k])/r.d[k],b=(hi[k]-r.o[k])/r.d[k];if(a>b)std::swap(a,b);near=std::max(near,a);far=std::min(far,b);if(far<near)return false;}return true;}};
bool intersect(Tri t,Ray r,double &distance,V &normal){V e=t.b-t.a,f=t.c-t.a,p=cross(r.d,f);double det=dot(e,p);if(std::abs(det)<1e-10)return false;V s=r.o-t.a;double u=dot(s,p)/det;if(u<0||u>1)return false;V q=cross(s,e);double v=dot(r.d,q)/det;if(v<0||u+v>1)return false;double d=dot(f,q)/det;if(d<=1e-6||d>=distance)return false;distance=d;normal=norm(cross(e,f));if(dot(normal,r.d)>0)normal=normal*(-1);return true;}
struct Node{Box box;int begin=0,end=0,left=-1,right=-1;};
#include "fast_bvh.h"
struct Scene{
 std::unique_ptr<FastBVH> fast;
 std::vector<Tri> tris;std::vector<Node> nodes;
#ifdef SOFTFRAME_EMBREE
 RTCDevice device=nullptr;RTCScene scene=nullptr;
 ~Scene(){if(scene)rtcReleaseScene(scene);if(device)rtcReleaseDevice(device);}
 void accelerate(){
  device=rtcNewDevice(nullptr);assert(device);scene=rtcNewScene(device);
  RTCGeometry geom=rtcNewGeometry(device,RTC_GEOMETRY_TYPE_TRIANGLE);
  struct Vertex{float x,y,z;};struct Index{unsigned a,b,c;};
  auto*v=(Vertex*)rtcSetNewGeometryBuffer(geom,RTC_BUFFER_TYPE_VERTEX,0,RTC_FORMAT_FLOAT3,sizeof(Vertex),tris.size()*3);
  auto*f=(Index*)rtcSetNewGeometryBuffer(geom,RTC_BUFFER_TYPE_INDEX,0,RTC_FORMAT_UINT3,sizeof(Index),tris.size());
  for(size_t i=0;i<tris.size();i++){auto t=tris[i];V xyz[3]={t.a,t.b,t.c};for(int j=0;j<3;j++)v[3*i+j]={float(xyz[j].x),float(xyz[j].y),float(xyz[j].z)};f[i]={unsigned(3*i),unsigned(3*i+1),unsigned(3*i+2)};}
  rtcCommitGeometry(geom);rtcAttachGeometry(scene,geom);rtcReleaseGeometry(geom);rtcCommitScene(scene);assert(rtcGetDeviceError(device)==RTC_ERROR_NONE);
 }
#else
 void accelerate(){fast=std::make_unique<FastBVH>(tris,nodes);}
#endif
 int build(int a,int b){
  int id=nodes.size();nodes.push_back({});Box box,cent;
  for(int i=a;i<b;i++){box.add(tris[i].a);box.add(tris[i].b);box.add(tris[i].c);cent.add(tris[i].center());}
  nodes[id].box=box;nodes[id].begin=a;nodes[id].end=b;if(b-a<=4)return id;
  auto area=[](const Box&b){V e=b.hi-b.lo;return 2*(e.x*e.y+e.x*e.z+e.y*e.z);};
  double best=1e300;int bestaxis=-1,bestsplit=-1;constexpr int bins=12;
  for(int axis=0;axis<3;axis++){
   double extent=cent.hi[axis]-cent.lo[axis];if(extent<1e-12)continue;
   Box boxes[bins];int counts[bins]{};
   for(int i=a;i<b;i++){int k=std::min(bins-1,int((tris[i].center()[axis]-cent.lo[axis])/extent*bins));counts[k]++;boxes[k].add(tris[i].a);boxes[k].add(tris[i].b);boxes[k].add(tris[i].c);}
   for(int split=0;split<bins-1;split++){Box left,right;int nl=0,nr=0;
    for(int k=0;k<bins;k++){if(!counts[k])continue;Box& dst=k<=split?left:right;dst.add(boxes[k].lo);dst.add(boxes[k].hi);(k<=split?nl:nr)+=counts[k];}
    if(!nl||!nr)continue;
    double cost=area(left)*nl+area(right)*nr;if(cost<best){best=cost;bestaxis=axis;bestsplit=split;}
   }
  }
  int m=(a+b)/2;
  if(bestaxis>=0){double scale=bins/(cent.hi[bestaxis]-cent.lo[bestaxis]);auto mid=std::partition(tris.begin()+a,tris.begin()+b,[&](const Tri&t){int k=std::min(bins-1,int((t.center()[bestaxis]-cent.lo[bestaxis])*scale));return k<=bestsplit;});m=mid-tris.begin();}
  if(m==a||m==b||bestaxis<0){m=(a+b)/2;std::nth_element(tris.begin()+a,tris.begin()+m,tris.begin()+b,[](const Tri&x,const Tri&y){return x.center().x<y.center().x;});}
  int l=build(a,m),r=build(m,b);nodes[id].left=l;nodes[id].right=r;return id;
 }

 bool trace(Ray ray,double &t,V &n,bool brute=false)const{bool found=false;
 if(!brute&&fast)return fast->trace(ray,t,n);
#ifdef SOFTFRAME_EMBREE
 if(!brute&&scene){RTCRayHit h{};h.ray.org_x=ray.o.x;h.ray.org_y=ray.o.y;h.ray.org_z=ray.o.z;h.ray.dir_x=ray.d.x;h.ray.dir_y=ray.d.y;h.ray.dir_z=ray.d.z;h.ray.tnear=1e-6;h.ray.tfar=t;h.ray.mask=~0u;h.hit.geomID=RTC_INVALID_GEOMETRY_ID;RTCIntersectArguments args;rtcInitIntersectArguments(&args);rtcIntersect1(scene,&h,&args);if(h.hit.geomID==RTC_INVALID_GEOMETRY_ID)return false;t=h.ray.tfar;n=norm({h.hit.Ng_x,h.hit.Ng_y,h.hit.Ng_z});if(dot(n,ray.d)>0)n=n*(-1);return true;}
#endif
 if(brute){for(auto tr:tris)found=intersect(tr,ray,t,n)||found;return found;}std::array<int,128> stack;int count=0;stack[count++]=0;while(count){auto &node=nodes[stack[--count]];if(!node.box.hit(ray,t))continue;if(node.left<0){for(int i=node.begin;i<node.end;i++)found=intersect(tris[i],ray,t,n)||found;}else{assert(count+2<=128);stack[count++]=node.left;stack[count++]=node.right;}}return found;}
};
#include "sobol.h"
struct RNG{
 uint64_t s;bool quasi=false;unsigned index=0,sample_number=0,dimension=0;unsigned shifts[32]{};
 void setup_quasi(){quasi=true;uint64_t h=s;for(auto &shift:shifts){h+=0x9e3779b97f4a7c15ULL;uint64_t z=h;z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;shift=unsigned(z^(z>>31));}}
 void sample_index(unsigned k){sample_number=k;index=k^(k>>1);dimension=0;}
 double get(){if(quasi){unsigned dim=dimension++,bits=index,value=shifts[dim%32];assert(dim<32);if(sample_number<512){value^=sobol_points[sample_number][dim];return (double(value)+.5)*0x1.0p-32;}while(bits){unsigned b=__builtin_ctz(bits);value^=sobol_directions[dim][b];bits&=bits-1;}return (double(value)+.5)*0x1.0p-32;}s^=s>>12;s^=s<<25;s^=s>>27;return ((s*2685821657736338717ULL)>>11)*0x1.0p-53;}
};
V basis(V n,V local){V t=norm(cross(std::abs(n.z)<.9?V{0,0,1}:V{0,1,0},n));return t*local.x+cross(n,t)*local.y+n*local.z;}
double D(double nh,double alpha){double a2=alpha*alpha,q=nh*nh*(a2-1)+1;return a2/(PI*q*q);}
double G1(double c,double a){return c<=0?0:2*c/(c+std::sqrt(a*a+(1-a*a)*c*c));}
// Mixture importance sampling: Lambertian + GGX NDF. Balance heuristic PDFs.
V bsdf(V n,V wo,V wi,double &pdf){double ni=std::max(0.,dot(n,wi)),no=std::max(0.,dot(n,wo));if(ni<=0||no<=0){pdf=0;return{};}V h=norm(wo+wi);double nh=std::max(0.,dot(n,h)),vh=std::max(1e-10,dot(wo,h)),alpha=.3;double dist=D(nh,alpha),f=.04+.96*std::pow(1-vh,5);pdf=.5*ni/PI+.5*dist*nh/(4*vh);return V{.6,.25,.12}*(1/PI)+V{1,1,1}*(f*dist*G1(ni,alpha)*G1(no,alpha)/(4*ni*no));}
V sample(V n,V wo,RNG &rng){double u=rng.get(),v=rng.get(),phi=2*PI*v;if(rng.get()<.5){double r=std::sqrt(u);return basis(n,{r*std::cos(phi),r*std::sin(phi),std::sqrt(1-u)});}double alpha=.3,ct=std::sqrt((1-u)/(1+(alpha*alpha-1)*u)),st=std::sqrt(1-ct*ct);V h=basis(n,{st*std::cos(phi),st*std::sin(phi),ct});return h*(2*dot(wo,h))-wo;}
V environment(V d){double sky=.3+.7*std::max(0.,d.y);return V{.7,.8,1}*sky;}
V path(const Scene&s,Ray r,RNG&rng){V throughput{1,1,1},sum{};for(int bounce=0;bounce<6;bounce++){double t=1e30;V n;if(!s.trace(r,t,n)){sum=sum+throughput*environment(r.d);break;}V p=r.o+r.d*t,wi=sample(n,r.d*(-1),rng);double pdf;V f=bsdf(n,r.d*(-1),wi,pdf);double cos=dot(n,wi);if(pdf<=0||cos<=0)break;throughput=throughput*f*(cos/pdf);if(bounce>=2){double survive=std::min(.95,std::max({throughput.x,throughput.y,throughput.z}));if(rng.get()>survive)break;throughput=throughput*(1/survive);}r={p+n*1e-4,wi};}return sum;}
struct Image{std::vector<V> pixels;uint64_t samples=0;double seconds;};
Image render(const Scene&s,int w,int h,int spp,bool adaptive,uint64_t seed=1234567,bool quasi=false){
 Image im;im.pixels.resize(w*h);auto start=std::chrono::steady_clock::now();
 std::vector<double> weights(w*h,1);double total=double(w*h);
 // Independent pilot samples choose a fixed per-pixel budget. Pilot values
 // are discarded: the final sample mean has no data-dependent stopping bias.
 if(adaptive){
  total=0;
#pragma omp parallel for reduction(+:total) schedule(dynamic,16)
  for(int i=0;i<w*h;i++){
   RNG rng{uint64_t(i+1)*(seed^0x9e3779b97f4a7c15ULL)};V mean{},m2{};
   for(int k=0;k<4;k++){V d=norm({(2*(i%w+rng.get())/w-1)*double(w)/h,1-2*(i/w+rng.get())/h,1.8});V v=path(s,{{0,0,-3},d},rng),delta=v-mean;mean=mean+delta*(1./(k+1));m2=m2+delta*(v-mean);}
   weights[i]=std::sqrt(std::max(1e-6,(m2.x+m2.y+m2.z)/9));total+=weights[i];
  }
  im.samples=uint64_t(w)*h*4;
 }
uint64_t final_samples=0;
#pragma omp parallel for reduction(+:final_samples) schedule(dynamic,16)
 for(int i=0;i<w*h;i++){
  RNG rng{uint64_t(i+1)*seed};if(quasi)rng.setup_quasi();V sum{};
  int count=adaptive?std::clamp(int(spp*(.8+.2*weights[i]*w*h/total)),4,spp*4):spp;
  for(int k=0;k<count;k++){if(quasi)rng.sample_index(k);V d=norm({(2*(i%w+rng.get())/w-1)*double(w)/h,1-2*(i/w+rng.get())/h,1.8});sum=sum+path(s,{{0,0,-3},d},rng);}
  im.pixels[i]=sum*(1./count);
  final_samples+=count;
 }
 im.samples+=final_samples;
 im.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();return im;
}
double error(const Image&a,const Image&b){double sum=0;for(size_t i=0;i<a.pixels.size();i++){V d=a.pixels[i]-b.pixels[i];sum+=dot(d,d);}return std::sqrt(sum/(3*a.pixels.size()));}
int main(){
#ifdef SOFTFRAME_EMBREE
 const char* backend="Embree 4 (external BVH)";
#else
 const char* backend="custom SIMD four-wide BVH";
#endif
 Scene s;RNG rng{42};for(int i=0;i<50000;i++){V a{rng.get()*20-10,rng.get()*20-10,rng.get()*20};s.tris.push_back({a,a+V{.05+rng.get()*.15,0,0},a+V{0,.05+rng.get()*.15,.03}});}s.build(0,s.tris.size());s.accelerate();
 for(int i=0;i<300;i++){Ray r{{0,0,-3},norm({rng.get()*2-1,rng.get()*2-1,1})};double a=1e30,b=1e30;V na,nb;assert(s.trace(r,a,na)==s.trace(r,b,nb,true));assert(std::abs(a-b)<1e-4);}
 std::vector<Ray> ray_batch(1000000);
 for(int i=0;i<1000000;i++){RNG local{uint64_t(i+1)*987654321};ray_batch[i]={{0,0,-3},norm({2*(i%1000+local.get())/1000-1,2*(i/1000+local.get())/1000-1,1})};}
 uint64_t hits=0;std::vector<double> ray_times;
 for(int trial=0;trial<5;trial++){hits=0;auto start=std::chrono::steady_clock::now();
#pragma omp parallel for reduction(+:hits)
  for(int i=0;i<1000000;i++){double t=1e30;V n;hits+=s.trace(ray_batch[i],t,n);}
  ray_times.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
 }
 std::sort(ray_times.begin(),ray_times.end());double seconds=ray_times[2];

 auto reference=render(s,96,64,1024,false,987654321),fixed=render(s,96,64,64,false),adaptive=render(s,96,64,64,true,1234567,true);
 std::vector<double> fixed_times{fixed.seconds};for(int trial=0;trial<4;trial++)fixed_times.push_back(render(s,96,64,64,false).seconds);std::sort(fixed_times.begin(),fixed_times.end());fixed.seconds=fixed_times[2];
 std::ofstream curve("results/softframe_quality.json");curve<<"{\"reference_spp\":1024,\"timing_repetitions\":5,\"timing_statistic\":\"median\",\"pilot_spp\":4,\"independent_pilot\":true,\"adaptive_sampler\":\"digitally shifted Sobol\",\"fixed_sampler\":\"IID\",\"fixed_rmse\":"<<error(fixed,reference)<<",\"fixed_seconds\":"<<fixed.seconds<<",\"candidates\":[";
 bool first=true;for(int budget:{16,24,32,48,64,96,128}){auto candidate=render(s,96,64,budget,true,1234567,true);std::vector<double> times{candidate.seconds};for(int trial=0;trial<4;trial++)times.push_back(render(s,96,64,budget,true,1234567,true).seconds);std::sort(times.begin(),times.end());candidate.seconds=times[2];if(!first)curve<<",";first=false;curve<<"{\"budget\":"<<budget<<",\"samples_including_pilot\":"<<candidate.samples<<",\"seconds\":"<<candidate.seconds<<",\"rmse\":"<<error(candidate,reference)<<"}";}curve<<"]}\n";
 std::ofstream out("results/softframe.ppm",std::ios::binary);out<<"P6\n96 64\n255\n";for(V c:adaptive.pixels){for(double a:{c.x,c.y,c.z})out.put(char(std::clamp(int(255*std::pow(std::max(0.,a),1/2.2)),0,255)));}
 std::cout<<"{\"backend\":\""<<backend<<"\",\"threads\":"<<omp_get_max_threads()<<",\"independent_reference_seed\":true,\"triangles\":50000,\"bvh_reference_checks\":300,\"intersection_rays\":1000000,\"ray_timing_repetitions\":5,\"ray_workload\":\"coherent jittered primary rays\",\"ray_generation_timed\":false,\"rays_per_second\":"<<1000000/seconds<<",\"hits\":"<<hits<<",\"fixed_seconds\":"<<fixed.seconds<<",\"adaptive_seconds\":"<<adaptive.seconds<<",\"fixed_rmse\":"<<error(fixed,reference)<<",\"adaptive_rmse\":"<<error(adaptive,reference)<<",\"fixed_samples\":"<<fixed.samples<<",\"adaptive_samples\":"<<adaptive.samples<<",\"reference_spp\":1024}\n";
 return 0;
}
