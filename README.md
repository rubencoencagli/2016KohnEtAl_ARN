%%
% Code to generate Figure 3 of 
% Kohn, COen-Cagli, Kanitscheide, Pouget
% Annu. Rev. Neurosci. 2016. 39:237–56
% doi: 10.1146/annurev-neuro-070815-013851
%
% 
%
%
%
%
%

%% Figure 3C
%% Perils of extrapolating small N results to large N.
clear all

CVG = [.4];
for c=1:numel(CVG)
    colG = [.4 .4 .4] - CVG(c);
    FI_LHC=[]; % Lin Carandini Harris 2015: homogeneous tuning, independent Negative-Binomial with multiplicative global gain
    FI_LHCdiff=[]; % Lin Carandini Harris 2015, plus differential correlations
    FIshuf_LHC=[];
    FIshuf_LHCdiff=[];
    FIasympt = 100;

    NN= 9*2.^[1:9]%2.^[4:11];%
    for nn=NN%
        nn
        tic
        K1 = nn; % # neurons pref ori
        k1 = [-pi:(2*pi/K1):(pi*(K1-1)/K1)]'; % neurons pref orientation
        a=350; %%%**** chosen to match the information for CVG=0 in Lin et al 2015
                
        ampl=1;
        wid=15*pi/180;
        muG = .3;
        sigmaG = CVG(c) * muG;
        h = 600;
        muA = 0.055;
        F = 1.7;
        
        %%% traslation-invariant tuning curves and derivatives
        f1 = circGaussian(0*pi/180,k1,wid,ampl);
        f2 = circGaussian(6*pi/180,k1,wid,ampl);
        df = f2-f1;
        %%% tuning curves and derivatives
        af1 = (a.*f1)*muG + h*muA;
        af2 = (a.*f2)*muG + h*muA;
        adf = (a.*df)*muG;
        if 0 %%% Lin-Harris-Carandini covariance matrix
            SLHC1 = sigmaG^2*((a.*f1)*(a.*f1)')  + F*diag(af1);
            SLHC2 = sigmaG^2*((a.*f2)*(a.*f2)')  + F*diag(af2);
            SLHC = 0.5*(SLHC1+SLHC2);
        else %%% Lin-Harris-Carandini covariance matrix, second option
            f12 = circGaussian(3*pi/180,k1,wid,ampl);
            af12 = (a.*f12)*muG + h*muA;
            SLHC = sigmaG^2*((a.*f12)*(a.*f12)')  + F*diag(af12);
        end
        %%% LHC plus differential
        SLHCdiff = SLHC + adf*adf'/FIasympt;% + eye(nn)*200;%
        toc
        
        tic
        %%% Fisher information
        FI_LHC = [FI_LHC; adf'*(SLHC\adf)];
        FI_LHCdiff = [FI_LHCdiff; adf'*(SLHCdiff\adf)];
        FIshuf_LHC = [FIshuf_LHC; adf'*((SLHC.*eye(nn))\adf)];
        FIshuf_LHCdiff = [FIshuf_LHCdiff; adf'*((SLHCdiff.*eye(nn))\adf)];
        toc
    end

    figure(1);
    hold on; axis square
    plot(NN,FI_LHC,'-','Color',colG);
    plot(NN,FIshuf_LHC,'--','Color',colG)
    set(gca,'xscale','log','yscale','log')

    figure(2);
    hold on; axis square
    plot(NN,FI_LHCdiff,'-','Color',colG);
    plot(NN,FIshuf_LHCdiff,'--','Color',colG)
    set(gca,'xscale','log','yscale','log')
end



%% Figure 3A,B
%% examples of mis-estimation of FI

clear all
FIin=[];
FIE=[];
FILP=[];
FIEdiag = [];
FILPest=[];

AA = 40;
NN= [10 20 50 100 500 1000 2500];%
for nn=NN%round(10.^([1:.1:3.25]))
    for aa=1:AA % smooth FI by resampling gains
        nn
        tic
        K1 = nn; % # neurons pref ori
        k1 = [-pi:(2*pi/K1):(pi*(K1-1)/K1)]'; % neurons pref orientation
        g = 20; % response gain of V1
        alpha=0; %baseline
        beta=1/(pi/4)^2;
        sigma1 = (g*0.005); % internal additive noise std
        
        %%% random amplitudes
        a = g * getLognAmplitudes(K1, .25); % random amplitude factors
        a(a>100)=100;
        %%% traslation-invariant tuning curves and derivatives
        [f, df, d2f] = vonMises(ones(size(a)),k1,0,alpha,beta,1);
        %%% heterogeneous tuning curves and derivatives
        af = a.*f;
        adf = a.*df;
        %%% correlation matrix before p-step
        C=zeros(nn);
        for j=1:nn
            fij = vonMises(ones(size(a)),k1,k1(j),alpha,beta,1);
            C(:,j) =  fij;
        end
        %%% our covariance before p-step
        SL = C  .* (a*a') * sigma1^2 ;
        %%% our covariance after p-step
        SLP = SL + diag(af);% + eye(nn)*200;%
        clear SL
        %%% ecker covariance, after p-step
        c0=.1;  CE = C*c0;
        SE = (CE + eye(size(C))*(1-c0)) .* sqrt(af*af') ;
        SEdiag = diag(diag(SE));
        %%% estimate covariance with missing entries
        toc
        
        tic
        %%% Fisher information
        [~, indk] = min(abs(k1));
        FIin=[FIin; -d2f(indk)/sigma1^2];
        % %     FILP = [FILP; (adf'*pinv(SLP)*adf)];
        % %     FIE = [FIE; (adf'*pinv(SE)*adf)];
        FILP = [FILP; (adf'*(SLP\adf))];
        FIE = [FIE; (adf'*(SE\adf))];
        tmp = (SEdiag\adf);
        FIEdiag = [FIEdiag; (adf'*tmp)^2 / (tmp'*SE*tmp)];
        toc
    end
end

%%% perturb covariance as if measured experimentally
CLP = SLP./(diag(SLP).^.5*diag(SLP)'.^.5);
Cest=eye(size(CLP));
for i=1:20:K1-20
    tmp=[];
    for j=i:i+19
        tmp=[tmp; diag(CLP,j)];
    end
    ind = randi(numel(tmp),1,numel(tmp));
    tmpest = tmp(ind);
    for j=i:i+19
        Cest = Cest + diag(tmpest(((j-i)*K1-(j-i)*i-sum(1:(j-i-1)))+[1:K1-j]),j);
        %     Cest = Cest + diag(tmpest(1:K1-i),i) + ...
        %         diag(tmpest((K1-i)+[1:K1-i-1]),i+1) + ...
        %         diag(tmpest((2*K1-2*i-1)+[1:K1-i-2]),i+2) +...
        %         diag(tmpest((3*K1-3*i-3)+[1:K1-i-3]),i+3) +...
        %         diag(tmpest((4*K1-4*i-6)+[1:K1-i-4]),i+4);
    end
end
Cest = Cest + triu(Cest,1)';

Qest = zeros(size(C));
rmean = 2*(sum(sum(triu(Cest))))/(K1*(K1-1));
tmp1=sqrt((1-rmean)*(1-rmean+rmean*K1));
for i=1:K1
    Qest(i,i) = sqrt(rmean+2*(1-rmean-tmp1)/K1)*(1+tmp1)/(rmean*sqrt(K1));
    for j=i+1:K1
        Qest(i,j) = sqrt(Cest(i,j)+2*(1-Cest(i,j)-sqrt((1-Cest(i,j))*(1-Cest(i,j)+Cest(i,j)*K1)))/K1)/sqrt(K1);
        Qest(j,i) = Qest(i,j);
    end
end
SLPest = Qest*Qest';
SLPest = SLPest -diag(diag(SLPest)) + eye(size(SLPest));
SLPest = diag(diag(SLP).^.5)*SLPest*diag(diag(SLP).^.5);

for nn=NN%round(10.^([1:.1:3.25]))
    for aa=1:AA % smooth FI by resampling gains
        tic
        K1 = nn;
        ind = randperm(NN(end),K1);
        SLPestnn = SLPest(ind,ind);
        FILPest = [FILPest; (adf(ind)'*(SLPestnn\adf(ind)))];
        toc
    end
end


BB=numel(FIin);
FIin = nanmean(reshape(FIin,AA,BB/AA),1);
FILP = nanmean(reshape(FILP,AA,BB/AA),1);
FIE = nanmean(reshape(FIE,AA,BB/AA),1);
FIEdiag = nanmean(reshape(FIEdiag,AA,BB/AA),1);
FILPest = nanmean(reshape(FILPest,AA,BB/AA),1);

figure;
subplot(1,2,1); hold on; axis square
plot(NN,FILP,'-k');
plot(NN,FILPest,'--b')
set(gca,'xscale','lin','yscale','lin','TickDir','out','XTick',[0:500:5000],'YTick',[0:500:5000],'xlim',[10 10^3.5],'ylim',[0 2000]); 
subplot(1,2,2); hold on; axis square
plot(NN,FIE,'-k');
plot(NN,FIEdiag,'--b')
set(gca,'xscale','lin','yscale','lin','TickDir','out','XTick',[0:500:5000],'YTick',[0:1000:5000],'xlim',[10 10^3.5],'ylim',[0 5000]); 



