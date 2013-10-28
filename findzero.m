function varargout=findzero(t,f)

s=sign(f);
sd=diff(s);
is=find(abs(sd)>0);
tz=t(is);
varargout{1}=tz;
varargout{2}=is;
