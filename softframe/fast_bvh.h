// Project-owned four-wide SIMD BVH. No external traversal library.
#include <xmmintrin.h>
#include <memory>
struct FastBVH {
 struct alignas(16) Wide {
  float lx[4],ly[4],lz[4],hx[4],hy[4],hz[4];
  int child[4],count;
 };
 struct alignas(16) Leaf {
  float ax[4],ay[4],az[4],ex[4],ey[4],ez[4],fx[4],fy[4],fz[4];
  int count;
 };
 std::vector<Wide> wide;
 std::vector<Leaf> leaves;
 int root;
 FastBVH(const std::vector<Tri>&tri,const std::vector<Node>&tree){root=make(0,tri,tree);}
 int make(int b,const std::vector<Tri>&tri,const std::vector<Node>&tree){
  const auto &n=tree[b];
  if(n.left<0){Leaf leaf{};leaf.count=n.end-n.begin;assert(leaf.count<=4);
   for(int i=0;i<leaf.count;i++){auto t=tri[n.begin+i];V e=t.b-t.a,f=t.c-t.a;leaf.ax[i]=t.a.x;leaf.ay[i]=t.a.y;leaf.az[i]=t.a.z;leaf.ex[i]=e.x;leaf.ey[i]=e.y;leaf.ez[i]=e.z;leaf.fx[i]=f.x;leaf.fy[i]=f.y;leaf.fz[i]=f.z;}
   leaves.push_back(leaf);return -int(leaves.size());
  }
  std::vector<int> children{n.left,n.right};
  while(children.size()<4){int selected=-1,largest=0;for(size_t i=0;i<children.size();i++){auto&c=tree[children[i]];if(c.left>=0&&c.end-c.begin>largest){selected=i;largest=c.end-c.begin;}}if(selected<0)break;int id=children[selected];children.erase(children.begin()+selected);children.insert(children.begin()+selected,{tree[id].left,tree[id].right});}
  int id=wide.size();wide.emplace_back();Wide out{};out.count=children.size();
  for(int i=0;i<out.count;i++){auto box=tree[children[i]].box;out.lx[i]=box.lo.x;out.ly[i]=box.lo.y;out.lz[i]=box.lo.z;out.hx[i]=box.hi.x;out.hy[i]=box.hi.y;out.hz[i]=box.hi.z;out.child[i]=make(children[i],tri,tree);}
  wide[id]=out;return id;
 }
 bool trace(Ray ray,double&distance,V&normal)const{
  const __m128 ox=_mm_set1_ps(ray.o.x),oy=_mm_set1_ps(ray.o.y),oz=_mm_set1_ps(ray.o.z);
  const __m128 dx=_mm_set1_ps(ray.d.x),dy=_mm_set1_ps(ray.d.y),dz=_mm_set1_ps(ray.d.z);
  const auto inverse=[](double d){return float(std::abs(d)<1e-30?std::copysign(1e30,d):1/d);};
  const __m128 ix=_mm_set1_ps(inverse(ray.d.x)),iy=_mm_set1_ps(inverse(ray.d.y)),iz=_mm_set1_ps(inverse(ray.d.z));
  struct Item{int id;float near;};Item stack[128];int size=1;stack[0]={root,0};bool found=false;const Leaf*hitLeaf=nullptr;int hitLane=0;float closest=distance;
  while(size){Item item=stack[--size];if(item.near>closest)continue;
   if(item.id<0){const auto &l=leaves[-item.id-1];
    __m128 ex=_mm_load_ps(l.ex),ey=_mm_load_ps(l.ey),ez=_mm_load_ps(l.ez),fx=_mm_load_ps(l.fx),fy=_mm_load_ps(l.fy),fz=_mm_load_ps(l.fz);
    __m128 px=_mm_sub_ps(_mm_mul_ps(dy,fz),_mm_mul_ps(dz,fy)),py=_mm_sub_ps(_mm_mul_ps(dz,fx),_mm_mul_ps(dx,fz)),pz=_mm_sub_ps(_mm_mul_ps(dx,fy),_mm_mul_ps(dy,fx));
    __m128 det=_mm_add_ps(_mm_add_ps(_mm_mul_ps(ex,px),_mm_mul_ps(ey,py)),_mm_mul_ps(ez,pz));
    __m128 inv=_mm_div_ps(_mm_set1_ps(1),det);
    __m128 sx=_mm_sub_ps(ox,_mm_load_ps(l.ax)),sy=_mm_sub_ps(oy,_mm_load_ps(l.ay)),sz=_mm_sub_ps(oz,_mm_load_ps(l.az));
    __m128 u=_mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(sx,px),_mm_mul_ps(sy,py)),_mm_mul_ps(sz,pz)),inv);
    __m128 qx=_mm_sub_ps(_mm_mul_ps(sy,ez),_mm_mul_ps(sz,ey)),qy=_mm_sub_ps(_mm_mul_ps(sz,ex),_mm_mul_ps(sx,ez)),qz=_mm_sub_ps(_mm_mul_ps(sx,ey),_mm_mul_ps(sy,ex));
    __m128 v=_mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(dx,qx),_mm_mul_ps(dy,qy)),_mm_mul_ps(dz,qz)),inv);
    __m128 t=_mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(fx,qx),_mm_mul_ps(fy,qy)),_mm_mul_ps(fz,qz)),inv),zero=_mm_setzero_ps();
    __m128 absdet=_mm_andnot_ps(_mm_set1_ps(-0.f),det);
    __m128 valid=_mm_and_ps(_mm_cmpgt_ps(absdet,_mm_set1_ps(1e-10)),_mm_and_ps(_mm_cmpge_ps(u,zero),_mm_and_ps(_mm_cmpge_ps(v,zero),_mm_cmple_ps(_mm_add_ps(u,v),_mm_set1_ps(1)))));
    valid=_mm_and_ps(valid,_mm_and_ps(_mm_cmpgt_ps(t,_mm_set1_ps(1e-6)),_mm_cmplt_ps(t,_mm_set1_ps(closest))));
    int mask=_mm_movemask_ps(valid)&((1<<l.count)-1);alignas(16)float ts[4];_mm_store_ps(ts,t);
    while(mask){int i=__builtin_ctz(mask);mask&=mask-1;if(ts[i]>=closest)continue;closest=ts[i];hitLeaf=&l;hitLane=i;found=true;}
   }else{const auto &n=wide[item.id];
    __m128 ax=_mm_mul_ps(_mm_sub_ps(_mm_load_ps(n.lx),ox),ix),bx=_mm_mul_ps(_mm_sub_ps(_mm_load_ps(n.hx),ox),ix);
    __m128 ay=_mm_mul_ps(_mm_sub_ps(_mm_load_ps(n.ly),oy),iy),by=_mm_mul_ps(_mm_sub_ps(_mm_load_ps(n.hy),oy),iy);
    __m128 az=_mm_mul_ps(_mm_sub_ps(_mm_load_ps(n.lz),oz),iz),bz=_mm_mul_ps(_mm_sub_ps(_mm_load_ps(n.hz),oz),iz);
    __m128 near=_mm_max_ps(_mm_set1_ps(1e-6),_mm_max_ps(_mm_min_ps(ax,bx),_mm_max_ps(_mm_min_ps(ay,by),_mm_min_ps(az,bz))));
    __m128 far=_mm_min_ps(_mm_set1_ps(closest),_mm_min_ps(_mm_max_ps(ax,bx),_mm_min_ps(_mm_max_ps(ay,by),_mm_max_ps(az,bz))));
    int mask=_mm_movemask_ps(_mm_cmple_ps(near,far))&((1<<n.count)-1);alignas(16)float ds[4];_mm_store_ps(ds,near);Item hits[4];int count=0;
    while(mask){int i=__builtin_ctz(mask);mask&=mask-1;hits[count++]={n.child[i],ds[i]};}
    // Push farthest first so the nearest child reduces the distance bound.
    for(int i=1;i<count;i++){auto x=hits[i];int j=i;while(j&&hits[j-1].near<x.near){hits[j]=hits[j-1];j--;}hits[j]=x;}
    assert(size+count<=128);for(int i=0;i<count;i++)stack[size++]=hits[i];
   }
  }
  if(found){const auto&l=*hitLeaf;int i=hitLane;normal=norm(cross({l.ex[i],l.ey[i],l.ez[i]},{l.fx[i],l.fy[i],l.fz[i]}));distance=closest;if(dot(normal,ray.d)>0)normal=normal*(-1);}return found;
 }
};
