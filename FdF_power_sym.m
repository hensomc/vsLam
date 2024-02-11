function [F,dFx,dFy,dFxy,dFxx,dFyy]=FdF_power_sym(Nmax)
% : Generate F and its derivative in symbolic form
% F=[1 x y x^2 x*y y^2 ...] = terms in the Pascal triabgle that includes 
% combination of terms (1,x,x^2 ...x^Namx)
%                and   (1,y,y^2 ...y^Namx)

syms x y
xn=x;yn=y;

%F(x,y)
for i=1:length(xn)
  FXY=[1];
  %loop over N terms
  for N=1:Nmax
    x=xn(i);y=yn(i);
    J=0:N;I=N-J;
    xx=x.^I;
    yy=y.^J;
    FXY=[FXY xx.*yy];
  end
  FF(i,:)=FXY;
end
F=FF;

%dF/dx
for i=1:length(xn)
  FXY=[0];
  x=xn(i);y=yn(i);
  for N=1:Nmax
    J=0:N;
    I=N-J;
    xx=I.*x.^(I-1);
    yt=y;
    yy=yt.^J;
    FXY=[FXY xx.*yy];
  end
  dFx(i,:)=FXY;
end

%dF/dy
for i=1:length(xn)
  FXY=[0];
  x=xn(i);y=yn(i);
  for N=1:Nmax
    J=0:N;
    I=N-J;
    xx=x.^I;
    yy=J.*y.^(J-1);
    FXY=[FXY xx.*yy];
  end
  dFy(i,:)=FXY;
end

%dFxx
for i=1:length(xn)
  FXY=[0];
  x=xn(i);y=yn(i);
  for N=1:Nmax
    J=0:N;
    I=N-J;
    xx=I.*(I-1).*x.^(I-2);
    yy=y.^J;
    FXY=[FXY xx.*yy];
  end
  dFxx(i,:)=FXY;
end

%dFxy
for i=1:length(xn)
  FXY=[0];
  x=xn(i);y=yn(i);
  for N=1:Nmax
    J=0:N;
    I=N-J;
    xx=I.*x.^(I-1);
    yy=J.*y.^(J-1);
    FXY=[FXY xx.*yy];
  end
  dFxy(i,:)=FXY;
end

%dFyy
for i=1:length(xn)
  FXY=[0];
  x=xn(i);y=yn(i);
  for N=1:Nmax
    J=0:N;
    I=N-J;
    xx=x.^I;
    yy=J.*(J-1).*y.^(J-2);
    FXY=[FXY xx.*yy];
  end
dFyy(i,:)=FXY;
end


