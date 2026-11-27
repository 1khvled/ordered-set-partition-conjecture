// Small exact rational exploration only (NOT a proof artifact).
const B=(n,k)=>{if(k<0||k>n)return 0n;let r=1n;for(let i=1;i<=k;i++)r=r*BigInt(n-i+1)/BigInt(i);return r};
function rows(j,N){ // a[n][m] as BigInt
  const a=Array.from({length:N+1},()=>[]);
  a[0][0]=1n;
  for(let n=1;n<=N;n++){
    a[n][0]=0n;
    for(let m=0;m<=Math.floor(n/j);m++){
      let s=0n;
      for(let r=1;r<=n;r++){
        const c=B(n,r);
        if(r<j) s+=c*(a[n-r][m]||0n);
        else s+=c*(a[n-r][m-1]||0n);
      }
      a[n][m]=s;
    }
  }
  return a;
}
function bareiss(M){
  const n=M.length; if(n===0)return 1n;
  let A=M.map(r=>r.slice()), sign=1n, prev=1n;
  for(let k=0;k<n-1;k++){
    if(A[k][k]===0n){let p=k+1;while(p<n&&A[p][k]===0n)p++;if(p===n)return 0n;[A[k],A[p]]=[A[p],A[k]];sign=-sign;}
    for(let i=k+1;i<n;i++)for(let jj=k+1;jj<n;jj++)A[i][jj]=(A[i][jj]*A[k][k]-A[i][k]*A[k][jj])/prev;
    for(let i=k+1;i<n;i++)A[i][k]=0n;
    prev=A[k][k];
  }
  return sign*A[n-1][n-1];
}
const comb=(N,k)=>{const arr=Array.from({length:N+1},(_,i)=>i);const go=(a,k)=>{if(k===0)return[[]];if(a.length<k)return[];const[h,...t]=a;return[...go(t,k-1).map(x=>[h,...x]),...go(t,k)]};return go(arr,k)};
function check(j,N,sizes){
  const a=rows(j,N); let neg=0, tested=0, first=null;
  for(const s of sizes){
    const R=comb(N,s);
    for(const I of R)for(const J of R){
      const M=I.map(n=>J.map(m=>a[n][m]||0n));
      const d=bareiss(M); tested++;
      if(d<0n){neg++;if(!first)first={I,J,d:d.toString()};}
    }
  }
  return {j,N,sizes,tested,neg,first};
}
const j=Number(process.argv[2]||2), N=Number(process.argv[3]||12), s=(process.argv[4]||'2,3').split(',').map(Number);
console.log(JSON.stringify(check(j,N,s)));
