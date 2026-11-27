// Small exact rational exploration only. Targeted minor search.
const B=(n,k)=>{if(k<0||k>n)return 0n;let r=1n;for(let i=1;i<=k;i++)r=r*BigInt(n-i+1)/BigInt(i);return r};
function rows(j,N){const a=Array.from({length:N+1},()=>[]);a[0][0]=1n;for(let n=1;n<=N;n++){a[n][0]=0n;for(let m=0;m<=Math.floor(n/j);m++){let s=0n;for(let r=1;r<=n;r++){const c=B(n,r);s+=(r<j?c*(a[n-r][m]||0n):c*(a[n-r][m-1]||0n));}a[n][m]=s;}}return a;}
function det(M){const n=M.length;if(n===0)return 1n;let A=M.map(r=>r.slice()),sg=1n,pv=1n;for(let k=0;k<n-1;k++){if(A[k][k]===0n){let q=k+1;while(q<n&&A[q][k]===0n)q++;if(q===n)return 0n;[A[k],A[q]]=[A[q],A[k]];sg=-sg;}for(let i=k+1;i<n;i++)for(let t=k+1;t<n;t++)A[i][t]=(A[i][t]*A[k][k]-A[i][k]*A[k][t])/pv;for(let i=k+1;i<n;i++)A[i][k]=0n;pv=A[k][k];}return sg*A[n-1][n-1];}
const comb=(N,k)=>{const go=(a,k)=>{if(k===0)return[[]];if(a.length<k)return[];const[h,...t]=a;return[...go(t,k-1).map(x=>[h,...x]),...go(t,k)]};return go(Array.from({length:N+1},(_,i)=>i),k)};
function run(j,N){
  const a=rows(j,N);const R2=comb(N,2),R3=comb(N,3);
  let neg2=0,neg3=0,f=null;
  for(const I of R2)for(const J of R2){const d=det(I.map(n=>J.map(m=>a[n][m]||0n)));if(d<0n){neg2++;if(!f)f={I,J,d:d.toString()};}}
  for(const I of R3)for(const J of R3){const d=det(I.map(n=>J.map(m=>a[n][m]||0n)));if(d<0n){neg3++;if(!f)f={I,J,d:d.toString()};}}
  // random size 4,5
  let neg4=0,neg5=0;
  const rnd=k=>{const s=new Set();while(s.size<k)s.add(Math.floor(Math.random()*(N+1)));return [...s].sort((x,y)=>x-y)};
  for(let it=0;it<30000;it++){const I=rnd(4),J=rnd(4);if(det(I.map(n=>J.map(m=>a[n][m]||0n)))<0n){neg4++;if(!f)f={I,J,size:4};}}
  for(let it=0;it<10000;it++){const I=rnd(5),J=rnd(5);if(det(I.map(n=>J.map(m=>a[n][m]||0n)))<0n){neg5++;if(!f)f={I,J,size:5};}}
  console.log(JSON.stringify({j,N,neg2,neg3,neg4,neg5,first:f}));
}
run(Number(process.argv[2]||2),Number(process.argv[3]||18));
