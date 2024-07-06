%functions to include in the folder:
% condensate_substrate
% solver_FDCH
% solverParser
% condensate_substrate_pp (done)
% droppullet_package (done)
% Flory_Huggins_setup (done)
% phaseEqm (done)
% reaction (done)
% makegrid (done)


%A wrapper function to conduct phase field simulations in the VECTOR paper
run_base_case = false;
run_study = true; %perform simulations
pp = false; %do post-processing
movie_time = true; %perform simulations at time resolution needed for exporting movies
movie_save = false; %save movie
show_traj = false; %show substrate trajectory plots

params = [];
%diffusivities
Dab = 0.01;
Da = 10*Dab;
Db = 10*Dab;
Dc = Dab;
%chi parameters in Flory Huggins
chi_ab = 3;
chi_c = 3;
chi_ab_c = 2;
%size parameter in Flory Huggins
ra = 1;
rb = 1;
rab = ra+rb;
r = [ra,rb,rab];
rc = 2;
dof = 3;
%A: Core, B: IDR, AB: Corelet, C: telomere binding protein, S: Solvent
params.species = {'A','B','AB','C','S'};
params.metaspecies.At = [1,0,ra/rab,0];
params.metaspecies.Bt = [0,1,rb/rab,0];
params.chi = {'AB','S',chi_ab,'C','S',chi_c,'AB','C',chi_ab_c};
params.r = [ra,rb,rab,rc,1];
params.D0 = diag([Da,Db,Dab,Dc]);
params = Flory_Huggins_setup(params);
lambda = 0.03;
params.kappa = lambda^2;

Rcore = 0.273; %radius of the Corelet droplet R0, chosen to be the characteristic length to define characteristic time scals. Rcore/lambda = 9
td = Rcore^3/(Dab*sqrt(params.kappa)); %diffusional time scale

tr_to_td = 2000; %Da number for base case simulation
%base case dissolution k
k_base = Rcore^2/params.kappa/(td*tr_to_td);
%basis case friction
td_to_ts = 50;
friction_base = td/td_to_ts/Rcore*sqrt(params.kappa);

Kd_bright = 0;
Kd_dark = 0.2;
k_dark = k_base/Kd_dark;
k_bright_to_k_dark = 4e3;
k_bright = k_dark*k_bright_to_k_dark;

L = 1.1; %L/lambda = 36
%r1: radius of corelet when there are two corelet
%r2: radius of corelet when two merge into one
%r0: radius of telomere
%L/2 = 2r1 + 2*(2*r0)
%r2 = sqrt(2) * r1, r2=3*r0
r2_to_r0 = 2;
r0 = L/2 / (2/sqrt(2)*r2_to_r0 + 4);
r2 = r2_to_r0*r0;
r1 = r2/sqrt(2);

%get an estimate of the binodal when A and B bind to form AB, the Corelet
%three phase
%gamma: phase fraction
%gamma1(1): fraction of AB (corelet) droplet
bp.species = {'AB','C','S'};
bp.chi = {'AB','S',chi_ab,'C','S',chi_c,'AB','C',chi_ab_c};
bp.r = [rab,rc,1];
bp = Flory_Huggins_setup(bp);
psi = [0.8,0.1,0.1; 0.1,0.8,0.1];
[phi1,gamma,psi] = phaseEqm(bp,[0.3,0.3],psi);
gamma1 = [2*pi*r1^2/L^2 ; 2*pi*r0^2/L^2];
%the following set of parameters work well
gamma1(1) = gamma1(1)*2.5; %increase fraction of gamma1
%initial center of telomere
x0 = 0.40;
%this set of parameters result in (merging time)/(before merging time) = 50/200 (close to ratio in experiment), and after merging corelet radius is Rcore = 0.273, resulting in Rcore/lambda=9.1
%Note that this result is dependent on td so let's fix Dab.
gamma1 = [gamma1; 1-sum(gamma1)];
phiavg = sum(phi1.*gamma1,1);
%get an estimate of the binodal when A and B do not bind (Corelet does not yet form). So there is only two phases: telomere binding protein and the dilute phase
bp.species = {'A','C','S'};
bp.chi = {'C','S',chi_c};
bp.r = [ra,rc,1];
bp = Flory_Huggins_setup(bp);
psi = [0.1,0.95,0.1];
[phi2,gamma2,psi] = phaseEqm(bp,phiavg(1:end-1),psi);


%initial coordinates of telomere loci
tel = [-x0,0; x0,0];
params.L = L*[1,1];
params.N = 128*[1,1];
% params.method = 'diag'; %not used because copyright issue does not support the use of Krylov method
params.message = true;
[x,y] = makegrid(params.L,params.N);

if movie_time
  Nt = 300;
else
  Nt = 10;
end

%the following studies are available: 'base_case','friction_sym','friction_asym_1','friction_asym_2','stiffness_asym_1','stiffness_asym_2','KV_asym_1','KV_asym_2','Jeffreys_asym_1','Jeffreys_asym_2','Jeffreys_asym_mainfig','KV_Jeffreys_recoil','Rouse_sym'
%'base case' consists of three stages: 1. equilibration from the initial condition when the light is off, 2. equilibriation with the light on, until the Corelet droplet forms and merges into one, 3. turn off the light.
%the end of the second stage of base case will be used as initial condition for the other studies
%Therefore, if the initial condition is not available must run base case first.

study = 'Jeffreys_asym_mainfig';

if run_base_case
  study = 'base_case';
end

%tr: reaction time scale, td: diffusional time scale, ts: characteristic time scale of viscously-dominated motion of telomere motion (note that it is denoted by t_v in the paper)
%tr_to_td means tr/td, etc
%fl: dashpot's friction coefficient of left telomere, fr: dashpot's friction coefficient of right telomere  
%El: spring's dimensionless stiffness parameter of left telomere, Er: spring's dimensionless stiffness parameter of right telomere
%sl: spring's stiffness parameter of left telomere, sr: spring's stiffness parameter of right telomere
switch study
case 'friction_sym'
  %symmetric Newtonian model, study how dissolution time and capillary force changes with ts/td
  tr_to_td = 2000;
  ts_to_td = 0.02:0.02:0.1;
  fl = td*ts_to_td/Rcore*sqrt(params.kappa);
  fr = fl;
  tspan = linspace(0,200,Nt);
case 'friction_asym_1'
  %asymmetric Newtonian model, fix the friction of the left telomere and vary the right telomere
  tr_to_td = 10000;
  tsl_to_td = 0.02*ones(1,5);
  tsr_to_td = 0.02*(50:50:250);
  fl = td*tsl_to_td/Rcore*sqrt(params.kappa);
  fr = td*tsr_to_td/Rcore*sqrt(params.kappa);
  tspan = linspace(0,500,Nt);
case 'friction_asym_2'
  %asymmetric Newtonian model, fix the ratio of the friction of the left and right telomere and vary both
  tr_to_td = 10000;
  tsl_to_td = 0.03*(2.^(0:4));
  tsr_to_td = tsl_to_td*5;
  fl = td*tsl_to_td/Rcore*sqrt(params.kappa);
  fr = td*tsr_to_td/Rcore*sqrt(params.kappa);
  tspan = linspace(0,600,Nt);
case 'stiffness_asym_1'
  %asymmetric Kevin-Voigt model, fix the dashpot's friction coefficient, fix the spring stiffness of the left telomere and vary the right telomere
  tr_to_td = 10000;
  %nondimensionalized stiffness parameter E = R0*k/sigma
  El = 0.05*ones(1,5);
  Er = El.*(2:2:10);
  sl = El/(Rcore/sqrt(params.kappa));
  sr = Er/(Rcore/sqrt(params.kappa));
  ts_to_td = 1e-4;
  %friction in Kevin-Voigt model
  ff = td*ts_to_td/Rcore*sqrt(params.kappa);
  tspan = linspace(0,600,Nt);
case 'stiffness_asym_2'
  %asymmetric Kevin-Voigt model, fix the dashpot's friction coefficient, fix the ratio of the spring stiffness of the left and right telomere and vary both
  tr_to_td = 10000;
  %nondimensionalized E = R0*k/sigma
  El = 0.01*(2.^(0:4));
  Er = El*5;
  sl = El/(Rcore/sqrt(params.kappa));
  sr = Er/(Rcore/sqrt(params.kappa));
  ts_to_td = 1e-4;
  %friction in Kevin-Voigt model
  ff = td*ts_to_td/Rcore*sqrt(params.kappa);
  tspan = linspace(0,600,Nt);
case 'KV_asym_1'
  %asymmetric Kevin-Voigt model, fix the retardation time tau (tau/td = ts/td /E = 2), fix the spring stiffness of the left telomere and vary the right telomere
  tr_to_td = 10000;
  tsl_to_td = 0.1*ones(1,5);
  tsr_to_td = tsl_to_td.*(2:2:10);
  %nondimensionalized E = R0*k/sigma
  El = 0.05*ones(1,5);
  Er = El.*(2:2:10);
  sl = El/(Rcore/sqrt(params.kappa));
  sr = Er/(Rcore/sqrt(params.kappa));
  %friction in Kevin-Voigt model
  fl = td*tsl_to_td/Rcore*sqrt(params.kappa);
  fr = td*tsr_to_td/Rcore*sqrt(params.kappa);
  tspan = linspace(0,600,Nt);
case 'KV_asym_2'
  %asymmetric Kevin-Voigt model, fix the retardation time tau (tau/td = ts/td /E = 2), fix the ratio of the spring stiffness of the left and right telomere and vary both
  tr_to_td = 10000;
  %nondimensionalized E = R0*k/sigma
  El = 0.008*(2.^(0:4));
  Er = El*5;
  tsl_to_td = 2*El; %so tau/td = t_s/td /E = 2
  tsr_to_td = 2*Er;
  sl = El/(Rcore/sqrt(params.kappa));
  sr = Er/(Rcore/sqrt(params.kappa));
  %friction in Kevin-Voigt model
  fl = td*tsl_to_td/Rcore*sqrt(params.kappa);
  fr = td*tsr_to_td/Rcore*sqrt(params.kappa);
  tspan = linspace(0,600,Nt);
case 'Jeffreys_asym_1'
  %asymmetric Jeffreys model, fix the retardation time tau (tau/td = ts/td /E = 2), fix the Maxwell relaxation time, fix the spring stiffness of the left telomere and vary the right telomere
  tr_to_td = 10000;
  tsl_to_td = 0.1*ones(1,5);
  tsr_to_td = tsl_to_td.*(2:2:10);
  %nondimensionalized E = R0*k/sigma
  El = 0.05*ones(1,5);
  Er = El.*(2:2:10);
  sl = El/(Rcore/sqrt(params.kappa));
  sr = Er/(Rcore/sqrt(params.kappa));
  %friction in Kevin-Voigt model
  fl = td*tsl_to_td/Rcore*sqrt(params.kappa);
  fr = td*tsr_to_td/Rcore*sqrt(params.kappa);
  %time constant of Maxwell relaxation time
  tau_to_td = 2;
  tau = td*2;
  tspan = linspace(0,600,Nt);
case 'Jeffreys_asym_2'
  %asymmetric Jeffreys model, fix the retardation time tau (tau/td = ts/td /E = 2), fix the Maxwell relaxation time,  fix the ratio of the spring stiffness of the left and right telomere and vary both
  tr_to_td = 10000;
  %nondimensionalized E = R0*k/sigma
  El = 0.008*(2.^(0:4));
  Er = El*5;
  tsl_to_td = 2*El; %so tau/td = t_s/td /E = 2
  tsr_to_td = 2*Er;
  sl = El/(Rcore/sqrt(params.kappa));
  sr = Er/(Rcore/sqrt(params.kappa));
  %friction in Kevin-Voigt model
  fl = td*tsl_to_td/Rcore*sqrt(params.kappa);
  fr = td*tsr_to_td/Rcore*sqrt(params.kappa);
  %time constant of Maxwell relaxation time (denoted by tau_M in the paper)
  tau_to_td = 2;
  tau = td*2;
  tspan = linspace(0,600,Nt);
case 'Jeffreys_asym_mainfig'
  %used for main figure
  %similar to Jeffreys_asym_1,
  tr_to_td = 10000;
  ratio = [1,3,5];
  tsl_to_td = 0.1*ones(1,length(ratio));
  tsr_to_td = tsl_to_td.*ratio;
  %nondimensionalized E = R0*k/sigma
  El = 0.05*ones(1,length(ratio));
  Er = El.*ratio;
  sl = El/(Rcore/sqrt(params.kappa));
  sr = Er/(Rcore/sqrt(params.kappa));
  %friction in Kevin-Voigt model
  fl = td*tsl_to_td/Rcore*sqrt(params.kappa);
  fr = td*tsr_to_td/Rcore*sqrt(params.kappa);
  %time constant of Maxwell relaxation time (denoted by tau_M in the paper)
  tau_to_td = 2;
  tau = td*2;
  tspan = linspace(0,700,Nt);
case 'KV_Jeffreys_recoil'
  %compare the recoil dynamics of KV and Jeffreys by changing the Maxwell relaxation time
  %symmetric
  tr_to_td = 10000;
  %previously E=0.5, ts_to_td=2*E
  E = 0.4;
  ts_to_td = 2*E; %so tau/td = t_s/td /E = 2
  ss = E/(Rcore/sqrt(params.kappa));
  ff = td*ts_to_td/Rcore*sqrt(params.kappa);
  %inverse time constant of Maxwell relaxation time
  inv_tau_to_td = 0:0.2:0.8;
  tau = td./inv_tau_to_td;
  Nt = 300;
  tspan = linspace(0,1200,Nt);
case 'Rouse_sym'
  %symmetric Rouse
  %vary t_rs / td, t_rs is the Rouse relaxation time
  trs_to_td = 10.^(-4:0);
  srke = sqrt(td*trs_to_td/(4*pi)) * sqrt(params.kappa)/Rcore; %square root of kappa times eta
  tprocess = 500; %a characteristic time scale of the process. used to define kappa_to_eta so that tprocess * kappa_to_eta = 1
  kappa_to_eta = 1/tprocess;
  ff = srke / sqrt(kappa_to_eta);
  Nt = 600;
  tspan = linspace(0,1200,Nt);
end
if ~run_base_case
  k = Rcore^2/params.kappa/(td*tr_to_td);
  Kd_bright = 0;
  Kd_dark = 0.2;
  k_dark = k/Kd_dark;
  k_bright_to_k_dark = 4e3;
  k_bright = k_dark*k_bright_to_k_dark;
end


if show_traj
  figyadd = figure; %plot substrate motion
  tiledlayout('flow');
end

telcolor = [1,0,1];
corecolor = [0,1,0];

filepath = 'data/';

switch study
case 'base_case'
  Ncase = 1;
case {'stiffness_asym_1','stiffness_asym_2','KV_asym_1','KV_asym_2','Jeffreys_asym_1','Jeffreys_asym_2','Jeffreys_asym_mainfig'}
  Ncase = length(sr);
case {'friction_sym','friction_asym_1','friction_asym_2'}
  Ncase = length(fl);
case {'KV_Jeffreys_recoil'}
  Ncase = length(tau);
case 'Rouse_sym'
  Ncase = length(ff);
end

% the coordinates of the telomere over time
yadd = [];
% total force (capillary and viscoelastic) on the telomere over time
force = [];
% capillary force on the telomere over time
force_capillary = [];

for i = 1:Ncase

  if ~run_base_case
    vl = load([filepath,'initial_condition']);
    y0 = vl.y0;
  end

  params.substrate.dynamics = 'dashpot';
  params.substrate.U0 = [0,0,0,-1]; %telomere binding protein is attracted to the substrate
  params.substrate.kernel.func = 'Gaussian';
  params.substrate.kernel.params = sqrt(params.kappa*2);
  params.substrate.kernel.cutoff = 3*params.substrate.kernel.params;
  if run_base_case
    params.substrate.friction = friction_base;
  else
    switch study
    case {'friction_sym','friction_asym_1','friction_asym_2'}
      params.substrate.dynamics = 'dashpot';
      params.substrate.friction = [fl(i),fr(i)]; %asymmetric
    case {'stiffness_asym_1','stiffness_asym_2'}
      params.substrate.dynamics = 'Kevin-Voigt';
      params.substrate.friction = ff;
      params.substrate.spring.stiffness = [sl(i),sr(i)];
      %get origin from yadd after equilibration
      origin = y0(end-3:end);
      origin = reshape(origin,2,2);
      params.substrate.spring.origin = origin;
    case {'KV_asym_1','KV_asym_2'}
      params.substrate.dynamics = 'Kevin-Voigt';
      params.substrate.spring.stiffness = [sl(i),sr(i)];
      params.substrate.friction = [fl(i),fr(i)];
      %get origin from yadd after equilibration
      origin = y0(end-3:end);
      origin = reshape(origin,2,2);
      params.substrate.spring.origin = origin;
    case {'Jeffreys_asym_1','Jeffreys_asym_2','Jeffreys_asym_mainfig'}
      params.substrate.dynamics = 'SFM';
      params.substrate.spring.stiffness = [sl(i),sr(i)];
      params.substrate.friction = [fl(i),fr(i)];
      %get origin from yadd after equilibration
      origin = y0(end-3:end);
      origin = reshape(origin,2,2);
      params.substrate.spring.origin = origin;
      params.substrate.tau = tau;
    case {'KV_Jeffreys_recoil'}
      if isinf(tau(i))
        params.substrate.dynamics = 'Kevin-Voigt';
      else
        params.substrate.dynamics = 'SFM';
        params.substrate.tau = tau(i);
      end
      params.substrate.spring.stiffness = ss;
      params.substrate.friction = ff;
      origin = y0(end-3:end);
      origin = reshape(origin,2,2);
      params.substrate.spring.origin = origin;
    case {'Rouse_sym'}
      params.substrate.dynamics = 'Rouse';
      params.substrate.friction = ff(i);
      params.substrate.Rouse.kappa_to_eta = kappa_to_eta;
    end
  end
  params.substrate.location = tel;
  % params.substrate.method = 'diag'; %not used because copyright issue does not support the use of Krylov method
  params.grid = {x,y};
  params.origin = [min(x(:)),min(y(:))];
  params = condensate_substrate_direct(params);

  %use for equilibration in the dark
  Kd_eqm = Kd_dark*ones(params.N);
  k_eqm = k_dark*ones(params.N);

  params_eqm = params;
  params_eqm.R = @(mu,c) reaction(mu,c,Kd_eqm(:),k_eqm(:),r);


  switch study
  case 'friction_sym'
    filename = [filepath,study,'_ts_to_td_',num2str(ts_to_td(i))];
  case {'friction_asym_1','friction_asym_2'}
    filename = [filepath,study,'_tsl_to_td_',num2str(tsl_to_td(i)),'_tsr_to_td_',num2str(tsr_to_td(i))];
  case {'stiffness_asym_1','stiffness_asym_2'}
    filename = [filepath,study,'_El_',num2str(El(i)),'_Er_',num2str(Er(i))];
  case {'KV_asym_1','KV_asym_2'}
    filename = [filepath,study,'_El_',num2str(El(i)),'_Er_',num2str(Er(i)),'_tsl_to_td_',num2str(tsl_to_td(i)),'_tsr_to_td_',num2str(tsr_to_td(i))];
  case {'Jeffreys_asym_1','Jeffreys_asym_2','Jeffreys_asym_mainfig'}
    filename = [filepath,study,'_El_',num2str(El(i)),'_Er_',num2str(Er(i)),'_tsl_to_td_',num2str(tsl_to_td(i)),'_tsr_to_td_',num2str(tsr_to_td(i)),'_tau_to_td',num2str(tau_to_td)];
  case {'KV_Jeffreys_recoil'}
    filename = [filepath,study,'_E_',num2str(E),'_ts_to_td_',num2str(ts_to_td),'_td_to_tau',num2str(inv_tau_to_td(i))];
  case {'Rouse_sym','base_case'}
    filename = [filepath,study];
  end

  light = x.^2+y.^2<(x0)^2;
  Kd = zeros(size(x));
  Kd(~light) = Kd_dark;
  Kd(light) = Kd_bright;
  k = zeros(size(x));
  k(~light) = k_dark;
  k(light) = k_bright;
  Kd = Kd(:);
  k = k(:);
  params.R = @(mu,c) reaction(mu,c,Kd,k,r);

  g0 = generateDroplet(x,y,tel,lambda,gamma2(1));
  AandB0 = phi2(1,1)*g0 + phi2(2,1)*(1-g0);
  C0 = phi2(1,2)*g0 + phi2(2,2)*(1-g0);
  AB0 = 0.01*AandB0;
  A0 = (ra/rab) * (AandB0-AB0);
  B0 = (rb/rab) * (AandB0-AB0);

  %these are times for equilibration
  tspan_on = linspace(0,400,Nt);
  tspan_off = linspace(0,400*5/3,Nt*5/3);

  if run_base_case
    %equilibrate
    thisy0 = [A0(:);B0(:);AB0(:);C0(:)];
    nstage = 3;
    for k = 1:nstage
      switch k
      case 1
        params_eqb = params_eqm;
        thisy0 = [thisy0; tel(:)];
        tspan = tspan_on;
      case 2
        params_eqb = params;
        tspan = tspan_on;
      case 3
        params_eqb = params_eqm;
        tspan = tspan_off;
      end
      [t,yy] = solver_FDCH_direct(tspan,thisy0,params_eqb);
      thisy0 = yy(end,:).';
      if k==2
        %use the end of this stage as the initial condition for studies
        y0 = thisy0;
      end
    end
    tspan = [tspan_on,tspan_off+tspan_on(end)];
    yy = solverParser(cat(1,yraw{:}),params,{'At','C'},[corecolor;telcolor],'clim',[0.05,0.5;0.05,1]);
    if movie_save
      vid = VideoWriter(filename,'MPEG-4');
      vid.FrameRate = 30;
      open(vid);
      rgb = yy.rgb;
      rgb(rgb<0) = 0;
      rgb = permute(rgb,[2,1,3,4]);
      writeVideo(vid,rgb);
      close(vid);
    end
    return
  end


  %turn off light
  thisy0 = y0;
  if isequal(params.substrate.dynamics,'SFM')
    thisy0 = [thisy0(:); origin(:)];
  end
  if isequal(params.substrate.dynamics,'Rouse')
    tel_now = thisy0(end-3:end);
    thisy0 = thisy0(1:end-4);
    tel_now = reshape(tel_now,2,2);
    params_eqm.substrate.location = tel_now;
    params_eqm = condensate_substrate_direct(params_eqm);
    thisy0 = [thisy0; params_eqm.substrate.location(:)];
  end
  [t,yraw] = solver_FDCH_direct(tspan,thisy0,params_eqm);

  if isequal(params.substrate.dynamics,'Rouse')
    params_eqm.substrate.N = params_eqm.add.N/2;
  end
  yy = solverParser(yraw,params_eqm,{'At','C'},[corecolor;telcolor],'clim',[0.05,0.5;0.05,1]);
  rgb = yy.rgb;
  rgb = permute(rgb,[2,1,3,4]);
  if ~movie_time
    figure;
    montage(rgb);
  end
  if movie_save
    vid = VideoWriter(filename,'MPEG-4');
    vid.FrameRate = 30;
    open(vid);
    rgb(rgb<0) = 0;
    writeVideo(vid,rgb);
    close(vid);
  end

  if isequal(params.substrate.dynamics,'SFM')
    thisyadd = reshape(yy.add,[],2,2,2);
    yadd_x0(:,:,i) = thisyadd(:,:,1,2);
    thisyadd = thisyadd(:,:,1,1);
  elseif isequal(params.substrate.dynamics,'Rouse')
    thisyadd = reshape(yy.add,[],1+params.substrate.Rouse.N,2,2);
    thisyadd = squeeze(thisyadd(:,1,:,1));
  end
  yadd(:,:,i) = thisyadd;

  if show_traj
    figure(figyadd);
    nexttile;
    plot(t,thisyadd);
  end

  thisforce = condensate_substrate_pp(params_eqm,yraw,'net');
  force(:,:,i) = thisforce(:,:,1);
  thisforce = condensate_substrate_pp(params_eqm,yraw,'capillary');
  force_capillary(:,:,i) = thisforce(:,:,1);


  %generate figures for visualization using VECTOR_package
  if ismember(study,{'Jeffreys_asym_mainfig','KV_asym_1','Jeffreys_asym_1','Rouse_sym'})
    resultfilepath = 'figures/';
    clim = [0.05,0.5;0.05,1];
    distance = diff(thisyadd,[],2);
    distance = abs(distance);
    %define the moment of first contact to be the point of highest positive second derivative
    d2d = diff(distance,2,1);
    [~,contact] = max(d2d,[],1); %just before the moment of contact
    d2 = diff(thisyadd([1,contact],2)); %displacement of the second telomere
    [~,before] = min(abs(thisyadd(:,2)-thisyadd(1,2)-d2/3)); %another points of interest before contact that is defined by the moment when the second telomere has travelled 30% of its displacement
    t_on = 0;
    thisparams = params;
    if isequal(params.substrate.dynamics,'Rouse')
      thisparams.substrate.N = thisparams.add.N/2;
    else
      thisparams.add.N = 4; %force solverParser to parse only the first 4 additional dof
    end
    snapshots = [1, before, contact-2, size(distance,1)];
    if isequal(study,'Rouse_sym')
      snapshots = [1, before, contact-1, size(distance,1)/2];
      VECTOR_package(tspan(1:end/2),yraw(1:end/2,:),snapshots,td,lambda,thisparams,clim,t_on,resultfilepath,[study,'_',num2str(i)],'new');
    else
      VECTOR_package(tspan,yraw,snapshots,td,lambda,thisparams,clim,t_on,resultfilepath,[study,'_',num2str(i)],'new');
    end
    close all
  end

end

switch study
case 'friction_sym'
  save([filepath,'yadd_',study],'t','yadd','force','force_capillary','fl','fr','ts_to_td');
case {'friction_asym_1','friction_asym_2'}
  save([filepath,'yadd_',study],'t','yadd','force','force_capillary','fl','fr','tsl_to_td','tsr_to_td');
case {'stiffness_asym_1','stiffness_asym_2'}
  save([filepath,'yadd_',study],'t','yadd','force','force_capillary','sl','sr','El','Er','ff','ts_to_td');
case {'KV_asym_1','KV_asym_2'}
  save([filepath,'yadd_',study],'t','yadd','force','force_capillary','sl','sr','El','Er','fl','fr','tsl_to_td','tsr_to_td');
case {'Jeffreys_asym_1','Jeffreys_asym_2','Jeffreys_asym_mainfig'}
  save([filepath,'yadd_',study],'t','yadd','force','force_capillary','sl','sr','El','Er','fl','fr','tsl_to_td','tsr_to_td','tau','tau_to_td');
case {'KV_Jeffreys_recoil'}
  save([filepath,'yadd_',study],'t','yadd','yadd_x0','force','force_capillary','ss','E','ff','ts_to_td','tau','inv_tau_to_td');
case {'Rouse_sym'}
  save([filepath,'yadd_',study],'t','yadd','force','force_capillary','ff','trs_to_td','kappa_to_eta');
end



function [R,dRdmu,dRdc] = reaction(mu,c,Kd,k,r)
  stoic = [-r(1),-r(2),r(3),0];
  sz = size(mu);
  ns = sz(end);
  mu = reshape(mu,[],ns);
  Rf = exp(mu(:,1)*r(1)+mu(:,2)*r(2));
  Rb = Kd.*exp(mu(:,3)*r(3));
  R = Rf-Rb;
  R = k.*R;
  R = R*stoic;
  R = reshape(R,sz);
  if nargout>1
    dRdmu = permute([Rf*r(1),Rf*r(2),-Rb*r(3),zeros(size(Rf))],[1,3,2]) .* stoic;
    dRdmu = k.*dRdmu;
    dRdmu = reshape(dRdmu,[sz,ns]);
    dRdc = [];
  end
end