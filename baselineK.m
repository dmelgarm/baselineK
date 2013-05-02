function [d v a varargout]=baselineK(t,ai,method,varargin)

% D.Melgar 05/2012 - 11/2012
% Correct strong motion accelerogram for baseline offset. Remember to do
% zero order corrections (pre-event mean) BEFORE running this program!!! 
% The nomenclature and variable names follow that of the Boore and Wang papers 
% (see REFERENCES) so you can use those as a manual to understand what's 
% going on here.
%
% USAGE
% [d v a]=baselineK(t,ai,'BoorePieceWise',t1,t2)
% [d v a]=baselineK(t,ai,'BooreQuad',t1)
% [d v a]=baselineK(t,ai,'Filter',Fc)
% [d v a t2fit Km]=baselineK(t,ai,'Static',t1,t2,K)
% [d v a t2fit Km]=baselineK(t,ai,'StaticSearch',t1,t2,tsearch,K)
% [d v a t1fit t2fit]=baselineK(t,ai,'StepFit',tstart,maxE)
% [d v a t1fit t2fit]=baselineK(t,ai,'StepFit',tstart,maxE,K)
%
% INPUTS
% t - Time vector
% ai - Acceleration time series
% t1 - Baseline correction start time
% t2 - Baseline correction intermediate time
% tstart - Start of motion, usually a p-wave pick
% Fc - Filter corner frequency
% tsearch - Time in seconds to be used for averaging in search
% maxE - A value from 0 to 1, proportion of energy accumulated to be used
%          for computation of static trend
% K - Value of static field constraint
% method - Baseline correction method
%    * 'BoorePieceWise' ~ Fit Piece-wise linear function as per Boore (1999,2001)
%    * 'BooreQuad' ~ Fit quadratic as per Boore et al. (2002)
%    * 'Filter' ~ High pass filter with corner frequency Fc
%    * 'Static' ~ Use Static field constraints i.e. from GPS and find the best correction 
%                 such that the mean of the last segment of the record is closest to it.
%    * 'StaticSearch' ~ Same as above but use  a grid search to find the best fit
%    * 'StepFit' ~ This will use the Wang et al (2011) method of fitting a
%                  step function to find the best values of t1 and t2 it
%                  can be used both with and without a static field
%                  constraints. This is the best method for automatic
%                  computation (it is not perfect though go read the paper!)
%
% OUTPUTS
% d - Baseline corrected displacement
% v - Baseline corrected velocity
% a - Baseline corrected acceleration
% t2fit - Initial correction time of best fit
% t2fit - Intermediate correction time of best fit
% Km - Misfit between static field constraint and actual achieved final
%      offset. It can be exactly zero although this is not always 
%      achievable
%
% CITING THIS CODE
%
% Please reference
% D.Melgar, Bock, Y., Snachez, D. & Crowell, B.W. (2013). On Robust and Automated Baseline
% Corrections for Strong Motion Seismology, J. Geophys. Res, in press.
%
% REFERENCES
% * Boore, D. M. (1999). Effect of baseline corrections on response spectra for
% two recordings of the 1999 Chi-Chi, Taiwan, earthquake, U.S. Geol. Surv. 
% Open-File Rept. 99-€“545, 37 pp.
% * Boore, D. M. (2001). Effect of baseline corrections on displacements and
% response spectra for several recordings of the 1999 Chi-Chi, Taiwan, earthquake,
% Bull. Seismol. Soc. Am. 91, 1199-€“1211.
% * Boore, D. M., C. D. Stephens, and W. B. Joyner (2002). Comments on baseline
% correction of digital strong-motion data: Examples from the 1999 Hector Mine,
% California, earthquake, Bull. Seismol. Soc. Am. 92, 1543-1560.
% * Wang, R., Schurr, B., Milkereit, C., Shao, Z. & Jin, M (2011). An Improved 
% Automatic Scheme for Empirical Baseline Correction of Digital Strong Motion 
% Records, Bull Seism. Soc. Am., 101(5), 2029-2044.


%Parse command line
if strcmp(method,'BoorePieceWise')
    t1=varargin{1};
    t2=varargin{2};
    method=1;
elseif strcmp(method,'BooreQuad')
    t1=varargin{1};
    method=2;
elseif strcmp(method,'Filter')
    Fc=varargin{1};
    method=3;
elseif strcmp(method,'Static')
    t1=varargin{1};
    t2=varargin{2};
    K=varargin{3};
    method=4;
elseif strcmp(method,'StaticSearch')
    t1=varargin{1};
    t2=varargin{2};
    dt=t(2)-t(1);
    tsearch=floor(varargin{3}/dt);
    K=varargin{4};
    method=5;
elseif strcmp(method,'StepFit')
    dtg=1; %Spacing for grid search
    t1=varargin{1};
    maxE=varargin{2};
    method=6;
    if nargin==6  %Include static field constraint
        K=varargin{3};
    end
else
    display('FATAL ERROR: Unknown correction method, exiting...')
    return
end

%Preliminary stuff
%ai must be a column vector
if size(ai,1)<size(ai,2)
    ai=ai';
end
%t must be a column vector
if size(t,1)<size(t,2)
    t=t';
end
%Clean up NaN
i=find(~isnan(t));
inan=find(isnan(t));
t=t(i);
ai=ai(i);
%Integrate to velocity to evaluate baseline corrections
vi=cumtrapz(t,ai);
di=cumtrapz(t,vi);
tf=max(t);

%Boore piecewise correction
if method==1
    %Evaluate af baseline af (t2 to end of record)
    iaf=find(t>=t2);
    G=[ones(size(vi(iaf))) t(iaf)];
    %Least squares for fit
    m=lsqlin(G,vi(iaf));
    iaf=find(t>=t2);
    v0=m(1);
    af=m(2);
    %Evaluate baseline am
    iam=find(t>=t1 & t<t2);
    vf=v0+af*t2;
    am=vf/(t2-t1);
    v0am=-am*t1;
    v1=v0am+am*t(iam);
    v2=v0+af*t(iaf);
    %Remove baselines and integrate records
    a=ai;
    a(iam)=a(iam)-am;
    a(iaf)=a(iaf)-af;
    a(inan)=NaN;
    t(inan)=NaN;
    v=cumtrapz(t,a);
    d=cumtrapz(t,v);
end

%Boore quadratic correction
if method==2
    %Evaluate af baseline af (t2 to end of record)
    i=find(t>=t1);
    G=[t(i).^2 t(i) ones(size(v(i)))];  %Inversion Kernel
    F=[t1.^2 t1 1]; %Constraints
    h=0;
    x=[G'*vi(i);h];
    G=[G'*G F';F 0];
    %Least squares for fit
    m=lsqlin(G,x);
    p=m(1);
    q=m(2);
    r=m(3);
    %Evaluate baseline am
    am=2*p*t(i)+q;
    %Remove baselines and integrate records
    a=ai;
    a(i)=a(i)-am;
    a(inan)=NaN;
    t(inan)=NaN;
    v=cumtrapz(t,a);
    d=cumtrapz(t,v);
end

%High pass filter
if method==3
    dt=t(2)-t(1);
    Fc=Fc*dt;
    [p,q]=butter(4,Fc,'high');
    a=filtfilt(p,q,ai);
    a(inan)=NaN;
    v=cumtrapz(t,a);
    d=cumtrapz(t,v);
end
%Static field constraints
if method==4
    %Evaluate af baseline af (t2 to end of record)
    iaf=find(t>=t2);
    G=[ones(size(vi(iaf))) t(iaf)];
    %Least squares for fit
    m=lsqlin(G,vi(iaf));
    v0=m(1);
    af=m(2);
    %Find time
    kfit=find(t<=169+0.001 & t>=169-0.001)
    t2fit=(af*tf^2+2*v0*tf-v0*t1+2*(K-d(kfit)))/(v0+af*t1);
    if t2fit<t1 || t2fit>t2  %Can't find exact fit
        display('Error: Baseline correction time is out of bounds.')
        display(['  Allowed interval is (' num2str(t1) ',' num2str(t2) ')']);
        display(['  tfit = ' num2str(t2fit)])
        display('  Cannot find exact solution! Iterating to find best fitting time...')
        t2fit=linspace(t1,t2,200);  %Do it iteratively and find best one
        for k=1:length(t2fit)
            [dfit vfit afit]=baselineK(t,ai,'BoorePieceWise',t1,t2fit(k)); %A little recursion never hurt anyone
            dend(k)=dfit(end)-K;
        end
        [mind i]=min(abs(dend));
        t2fit=t2fit(i);
        [d v a]=baselineK(t,ai,'BoorePieceWise',t1,t2fit);
        d(inan)=NaN;
        v(inan)=NaN;
        a(inan)=NaN;
        varargout{1}=t2fit; %Time for best fit
        varargout{2}=dend(i); %Misfit with regards to static field constraint
        return
    end    
    iaf=find(t>=t2fit); 
    %Evaluate baseline am
    iam=find(t>=t1 & t<t2fit);
    vf=v0+af*t2fit;
    am=vf/(t2fit-t1);
    %Remove baselines and integrate records
    a=ai;
    a(iam)=a(iam)-am;
    a(iaf)=a(iaf)-af;
    a(inan)=NaN;
    t(inan)=NaN;
    v=cumtrapz(t,a);
    d=cumtrapz(t,v);
    varargout{1}=t2fit;  %Time for best fit
    varargout{2}=0;   %Misfit with regards to static field constraint
end

if method==5
    tend=max(size(t));
    t2=t(end)-tsearch*dt;
    t2fit=linspace(t1,t2,200);  %Do it iteratively and find best one
    for k=1:length(t2fit)
        [dfit vfit afit]=baselineK(t,ai,'BoorePieceWise',t1,t2fit(k)); %A little recursion never hurt anyone
        dend(k)=mean(dfit(tend-tsearch:tend))-K;
    end
    [mind i]=min(abs(dend));
    t2fit=t2fit(i);
    [d v a]=baselineK(t,ai,'BoorePieceWise',t1,t2fit);
    d(inan)=NaN;
    v(inan)=NaN;
    a(inan)=NaN;
    varargout{1}=t2fit; %Time for best fit
    varargout{2}=dend(i); %Misfit with regards to static field constraint
    return
end

if method==6
    tic
    tend=max(size(t));
    %Find bounds for t2 grid search
    [foo itpga]=max(abs(ai));
    tpga=t(itpga);
    [tz iz]=findzero(t,di);
    tz=tz(end);
    iz=iz(end);
    t2start=max([tpga tz]);
    %Get energy of the signal and find the time when maxE % of energy has been released
    E=cumtrapz(t,abs(ai).^2);
    iE=find(E>=maxE*E(end));
    t2end=t(iE(1));
    if t2end<t2start   %Use last 10% of record
        trange=t(end)-t(1);
        t2end=t(end)-0.1*trange;
        iE=find(t>t2end);
    end
    if t2end<t2start   %Use last 5% of record
        trange=t(end)-t(1);
        t2end=t(end)-0.05*trange;
        iE=find(t>t2end);
    end
    if t2end<t2start   %Use last 5% of record
        d=di;
        v=vi;
        a=ai;
        varargout{1}=0;
        varargout{2}=0
        return
    end
    %Do fit for vf and af
    %Do fit of straight line segment from tmaxE to end of record
    G=[ones(size(vi(iE))) t(iE)];
    m=lsqlin(G,vi(iE));
    v0=m(1);
    af=m(2);
    %Find bounds for t1 grid search
    [dpgd id0]=max(abs(di(1:iz)));
    t1start=t(id0);
    %Get parameters for grid search
    T2=t2start:dtg:t2end;
    totaliter=0;
    nt=length(T2);
    for p=1:nt
        display(['Iterating through '  num2str(p) ' out of ' num2str(nt) ' possible values of t2.'])
        t2=T2(p);
        T1=t1start:dtg:T2(p);
        for q=1:length(T1);
            totaliter=totaliter+1;
            t1=T1(q);
            t1all(totaliter)=t1;
            t2all(totaliter)=t2;
            %Note fit is from tmaxE to tend but baseline is removed from t2 to tend
            iaf=find(t>=t2);
            %Evaluate baseline am
            iam=find(t>=t1 & t<t2);
            vf=v0+af*t2;
            am=vf/(t2-t1);
            v0am=-am*t1;
            v1=v0am+am*t(iam);
            v2=v0+af*t(iaf);
            %Remove baselines and integrate records
            a=ai;
            a(iam)=a(iam)-am;
            a(iaf)=a(iaf)-af;
            a(inan)=NaN;
            t(inan)=NaN;
            v=cumtrapz(t,a);
            d=cumtrapz(t,v);
            %Compute fit to step function
            options=optimset('TolX',1e-9,'MaxFunEvals',5000,'MaxIter',5000,'Display','off');
            tH=decimate(t,100);
            dH=decimate(d,100);
            if nargin==5
                H=@(p,tH)p(1)./(1+exp(-2*1*(tH-p(2))));
                [hparams,resnorm(totaliter)] = lsqcurvefit(H,[0 T1(end)],tH,dH,[],[],options);
                resnorm(totaliter)=sum(abs(dH-H(hparams,tH)));
            else
                H=@(p,tH)K./(1+exp(-2*1*(tH-p(1))));
                [hparams,resnorm(totaliter)] = lsqcurvefit(H,T1(end),tH,dH,[],[],options);
                resnorm(totaliter)=sum(abs(dH-H(hparams,tH)));
            end
        end
    end
    %Now get best fitting one
    [resmin imin]=min(resnorm);
    t1=t1all(imin);
    t2=t2all(imin);
    iaf=find(t>=t2);
    %Evaluate baseline am
    iam=find(t>=t1 & t<t2);
    vf=v0+af*t2;
    am=vf/(t2-t1);
    v0am=-am*t1;
    v1=v0am+am*t(iam);
    v2=v0+af*t(iaf);
    %Remove baselines and integrate records
    a=ai;
    a(iam)=a(iam)-am;
    a(iaf)=a(iaf)-af;
    a(inan)=NaN;
    t(inan)=NaN;
    v=cumtrapz(t,a);
    d=cumtrapz(t,v);
    varargout{1}=t1;
    varargout{2}=t2;
end

    
    
    
