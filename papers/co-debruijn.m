(* ::Package:: *)

jj[2]=2
kk[2]=6
ll[2]=2
bb[t_]:=(2^t-(-1)^t)/3
flg[m_]:=Block[{r=0,mm=m},While[mm>1,mm=Floor[mm/2];r++];r]
jj[n_]:=jj[n]=Block[{t=0,q=n,lq},While[Mod[q,2]==0,q=q/2;t++];lq=flg[q];
           If[lq==0,2^(n-1)+bb[t],
             If[Mod[lq+If[Floor[4q/2^lq]==5,1,0],2]==0,
                  2^(n-1)+bb[t+2],2^(n-1)-bb[t+2]]]]
kk[n_]:=kk[n]=(2^(n-1)-2)(2^n+2)+2jj[n]+2
ll[n_]:=ll[n]=If[Mod[jj[n],4]==1,(2^(n-1)+1)(jj[n]-1),(2^(n-1)+1)(2^n+jj[n]-3)]
f[n_,k_]:=If[k==0,0,
            If[k>=2^n,f[n,Mod[k,2^n]],
             If[n==1,Error,
              If[n==2,Floor[k/2],
               If[Mod[n,2]==1,finc[n,k],fdub[n,k]]]]]]
finc[n_,k_]:=If[k>jj[n-1] && k<=2^(n-1)+jj[n-1],Mod[1+sigf[n-1,k+2^(n-1)-jj[n]],2],
                 Mod[sigf[n-1,k],2]]
fdub[n_,k_]:=If[ll[n/2]<kk[n/2],
                 If[k<=ll[n/2]+2,c[n/2,k-1],
                   If[k<=kk[n/2]+3,c[n/2,k-3],c[n/2,k-4]]],
                 If[k<=kk[n/2]+1,c[n/2,k-1],
                   If[k<=ll[n/2]+3,c[n/2,k-2],c[n/2,k-4]]]]
c[n_,k_]:=If[Mod[k,2]==0,fplus[n,k/2],fminus[n,(k-1)/2]]
fplus[n_,k_]:=If[k==0,0,
               If[k>=2^n+2,fplus[n,Mod[k,2^n+2]],
                 If[k<=jj[n]+1,f[n,k-1],f[n,k-2]]]]
fminus[n_,k_]:=If[k>=2^n-2,fminus[n,Mod[k,2^n-2]],
                If[k<jj[n],f[n,k+1],f[n,k+2]]]
sigf[n_,k_]:=If[k==0,0,
               If[k>=2^n,sigf[n,Mod[k,2^n]],
                 If[Mod[n,2]==1,sigfinc[n,k],sigfdub[n,k]]]]
sigfdub[n_,k_]:=If[n==2,If[k==3,1,0],
                 If[ll[n/2]<kk[n/2],
                  If[k<=ll[n/2]+2,sigc[n/2,k-1],
                     If[k<=kk[n/2]+3,1-sigc[n/2,k-3],sigc[n/2,k-4]]],
                  If[k<=kk[n/2]+2,sigc[n/2,k-1],
                     If[k<=ll[n/2]+3,1-sigc[n/2,k-2],sigc[n/2,k-4]]]]]
sigc[n_,k_]:=Mod[sigfplus[n,Ceiling[k/2]]+sigfminus[n,Floor[k/2]],2]
sigfplus[n_,k_]:=If[k==0,0,
                  If[k>=2^n+2,Floor[k/(2^n+2)]+sigfplus[n,Mod[k,2^n+2]],
                    If[k<=jj[n]+1,sigf[n,k-1],1+sigf[n,k-2]]]]
sigfminus[n_,k_]:=If[k==0,0,
                  If[k>=2^n-2,Floor[k/(2^n-2)]+sigfminus[n,Mod[k,2^n-2]],
                    If[k<jj[n],sigf[n,k+1],1+sigf[n,k+2]]]]
sigfinc[n_,k_]:=If[k==0,0,
                   If[jj[n-1]<k && k<=2^(n-1)+jj[n-1],
                       If[n>3,1,0]+k+sigsigf[n-1,k+2^(n-1)-jj[n]],sigsigf[n-1,k]]]
sigsigf[n_,k_]:=If[k==0,0,
                 If[n==2,Floor[Mod[k,8]/4],
                  If[k>=2^n,
                    If[ll[n/2]<kk[n/2],Floor[k/2^n],0]+sigsigf[n,Mod[k,2^n]],
                   If[ll[n/2]<kk[n/2],
                    If[k<=ll[n/2]+2,sigsigc[n/2,k-1],
                     If[k<=kk[n/2]+3,1+k+sigsigc[n/2,k-3],sigsigc[n/2,k-4]]],
                    If[k<=kk[n/2]+1,sigsigc[n/2,k-1],
                     If[k<=ll[n/2]+3,1+k+sigsigc[n/2,k-2],sigsigc[n/2,k-4]]]]]]]
sigsigc[n_,k_]:=If[Mod[k,2]==0,sigfplus[n,k/2],sigfminus[n,(k-1)/2]]
alf[g_]:=Table[g[[2k-1]],{k,Length[g]/2}]
bet[g_]:=Table[g[[2k]],{k,Length[g]/2}]
pure[g_]:=Block[{s=Apply[Plus,g]},s==0 || s==Length[g]]
solve[j_,k_,n_]:=If[j>=k,k+(2^n-2)Mod[2^(n-3)(j-k),2^(n-1)+1],
                      j+(2^n+2)Mod[2^(n-3)(k-j),2^(n-1)-1]]
posplus[n_,k_]:=k+1+If[k>jj[n],1,0]
posminus[n_,k_]:=k-1-If[k>jj[n],1,0]
getpos[g_]:=If[Mod[Length[g],2]==1,getodd[g],geteven[g]]
geteven[{0,0}]=0
geteven[{0,1}]=1
geteven[{1,1}]=2
geteven[{1,0}]=3
geteven[g_]:=Block[{a,b,n,j,k,jx,kx,lx},
              a=alf[g];b=bet[g];
              j=getpos[a];k=getpos[b];n=Length[a];
              If[pure[a]&&pure[b],purepos[n,a[[1]],b[[1]]],
                If[pure[a],jx=posplus[n,j];kx=posminus[n,k];
                 If[Mod[jx+kx,2]==1,If[a[[1]]==0,jx--,jx++]];
                 lx=2solve[jx,kx,n],
                If[pure[b],jx=posminus[n,j];kx=posplus[n,k];
                  If[Mod[jx+kx,2]==0,If[b[[1]]==0,kx--,kx++]];
                  lx=Mod[-1+2solve[kx,jx+1,n],2^(2n)-4],
                 jx=posplus[n,j];kx=posminus[n,k];
                 If[Mod[jx+kx,2]==0,lx=2solve[jx,kx,n],
                   jx=posminus[n,j];kx=posplus[n,k];
                   lx=-1+2solve[kx,jx+1,n]]]];
                lx+1+If[lx>=kk[n],1,0]+If[lx>=ll[n],2,0]]]
purepos[n_,ax_,bx_]:=If[ax==0,
                      If[bx==0,0,ll[n]+1+If[kk[n]<ll[n],1,0]],
                      If[bx==0,ll[n]+2+If[kk[n]<ll[n],1,0],jj[2n]]]
del[n_]:=2^n-jj[n+1]
getodd[g_]:=Block[{b,n,j,k},
               n=Length[g]-1;
               beta=Table[Mod[g[[k]]+g[[k+1]],2],{k,n}];
               j=geteven[beta];
               k=j-If[j>=jj[n],del[n],0]+If[j==jj[n] && del[n]==1,2^n,0];
               If[f[n+1,k]==g[[1]],k,k+
                If[j<jj[n],2^n-del[n],If[j==jj[n],del[n],2^n+del[n]]]]]
test[n_]:=Block[{ff},ff=Table[f[n,k],{k,2^n+n}];
                 Do[If[getpos[Table[ff[[p]],{p,k,k+n-1}]]!=k,Print[k]],{k,2^n-1}]]

getpos[{1,1,1,1,1,1}]



