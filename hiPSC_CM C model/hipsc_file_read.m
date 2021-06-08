%hiPSC model output file columns

% 1. sim_time
% 2. Vm
% 3. ical
% 4. ipca
% 5. incx
% 6. ina
% 7. inak
% 8. ikto
% 9. ik1
%10. iks
%11. ikur
%12. ikr
%13. jserca
%14. jrel
%15. Ki
%16. Nai
%17. Cai
%18. ikach
%19. if
%20. icab
%21. inab 
%22. jca_diff
%23. jtr
%24. Ca_up
%25. Ca_rel
%26. Ca_sub

apdcalc = 1; % make it 1 to calculate APDs

%ivfile = 'hipsc_IV.txt';
%ivfile = 'hipsc_IV_ko8.txt';
%ivfile = 'hipsc_IV_cao6.txt';
%ivfile = 'hipsc_IV_ncx0.txt';
ivfile = 'hipsc_IV_ikach3protocol.txt';

%ivfile = 'hipsc_IV_if0.5.txt';
%ivfile = 'hipsc_IV_ikur0.txt';
%ivfile = 'hipsc_IV_ikach5.txt';
%ivfile = 'hipsc_IV_isoif.txt';

clear neodata;

neodata=textread(ivfile);
fprintf('IV file read: %s\n',ivfile);

figure
plot(neodata(:,1),neodata(:,2))
title('Vm')

figure
subplot(331); plot(neodata(:,1),neodata(:,6)); title('INa')
subplot(332); plot(neodata(:,1),neodata(:,3)); title('ICaL')
subplot(333); plot(neodata(:,1),neodata(:,9)); title('IK1')
subplot(334); plot(neodata(:,1),neodata(:,12)); title('IKr')
subplot(335); plot(neodata(:,1),neodata(:,10)); title('IKs')
subplot(336); plot(neodata(:,1),neodata(:,11)); title('IKur')
subplot(337); plot(neodata(:,1),neodata(:,8)); title('Ito')
subplot(338); plot(neodata(:,1),neodata(:,7)); title('INaK')
subplot(339); plot(neodata(:,1),neodata(:,5)); title('INaCa')

figure 
subplot(311); plot(neodata(:,1),neodata(:,16)); title('Nai')
subplot(312); plot(neodata(:,1),neodata(:,15)); title('Ki')
subplot(313); plot(neodata(:,1),neodata(:,17)); title('Cai')

figure 
subplot(311); plot(neodata(:,1),neodata(:,24)); title('Ca_up')
subplot(312); plot(neodata(:,1),neodata(:,25)); title('Ca_rel')
subplot(313); plot(neodata(:,1),neodata(:,26)); title('Ca_sub')

%% pacemaking
% subplot(221); plot(neodata(:,1),neodata(:,3)); title('ICaL'); xlim([5000 7000]);
% subplot(222); plot(neodata(:,1),neodata(:,25)); title('ICaT');xlim([5000 7000]);
% subplot(223); plot(neodata(:,1),neodata(:,5)); title('NCX');xlim([5000 7000]);
% subplot(224); plot(neodata(:,1),neodata(:,26)); title('If');xlim([5000 7000]);

%---------------------------------------------------------------
% compute dv/dt
tt=neodata(:,1);
vv=neodata(:,2);
lv=length(tt); 

j1=1;
plus1indx=0;
dv1=0; 
clear dv1max dt1max
dv1max=0;

for i=6:lv   
%    dv1=(pt1(i)-pt1(i-1))/dt;
    if (vv(i)>=vv(i-1))
        dv1=(vv(i)-vv(i-1))/(tt(i)-tt(i-1));

    if ((dv1>dv1max(j1))&&(dv1>50))
            plus1indx = 1; 
            dv1max(j1)=dv1;
            dt1max(j1)=i;%i;
    end
    else
        if plus1indx ==1
            plus1indx =0;
            j1=j1+1;
            dv1max(j1)=0;
        end
    end
        
end
if(apdcalc)
fprintf('dvdtmax\t\t APD50\t\t APD70\t\t APD90\n');    
for i=1:j1-1
    
    dvdtmax=dv1max(i);
    if i==1
            minV=min(vv(1:dt1max(i+1)));
    else
            minV=min(vv(dt1max(i)-50:dt1max(i)));
    end
    if i ==j1-1
        maxV=max(vv(dt1max(i):lv));
    else
        maxV=max(vv(dt1max(i):dt1max(i+1)));
    end
    Vamp=maxV-minV;
    %maxpp=dt1max(i)+ find(vv(dt1max(i):dt1max(i+1))==maxV);
    maxpp=dt1max(i)+ find(vv(dt1max(i):dt1max(i)+5000)==maxV);
   
    
    v30=Vamp*0.50;
    amp30=maxV-v30;
    if i==j1-1
            a30 = find(vv(maxpp(1,1):length(vv))<=amp30);
    else
            a30 = find(vv(maxpp(1,1):dt1max(i+1))<=amp30);
    end
    APD30=tt(a30(1,1)+maxpp)-tt(dt1max(i)); % w.r.t. dv/dt max
    
    v50=Vamp*0.70;
    amp50=maxV-v50;
    if i==j1-1
            a50 = find(vv(maxpp(1,1):length(vv))<=amp50);
    else
            a50 = find(vv(maxpp(1,1):dt1max(i+1))<=amp50);
    end
    APD50=tt(a50(1,1)+maxpp)-tt(dt1max(i)); % w.r.t. dv/dt max
    
    v80=Vamp*0.90;
    amp80=maxV-v80;
    if i==j1-1
            a80 = find(vv(maxpp(1,1):length(vv))<=amp80);
    else
            a80 = find(vv(maxpp(1,1):dt1max(i+1))<=amp80);
    end
    APD80=tt(a80(1,1)+maxpp)-tt(dt1max(i)); % w.r.t. dv/dt max
    
    fprintf('%f\t%f\t%f\t%f\n',dvdtmax, APD30, APD50, APD80);
end
end % apdcalc

