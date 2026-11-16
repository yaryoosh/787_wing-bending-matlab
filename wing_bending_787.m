
clear; clc; close all;

%% USER SETTINGS

b         = 60.12;     % wingspan (m)
fusW      = 5.77;      % fuselage width (m)
sweep_deg = 32;        % quarter-chord sweep (deg)
MTOW      = 227900;    % max takeoff weight (kg)

c_root        = 11.0;  % chord at side of body (m)
c_kink        = 6.7;   % chord at the trailing-edge kink (m)
c_tip         = 1.75;  % tip chord (m)
yKink_fromSOB = 9.0;   % kink position outboard of side of body (m)
tc_root = 0.15;  tc_kink = 0.11;  tc_tip = 0.09;   % thickness/chord

boxDepthFrac = 0.85;   % box depth as fraction of max wing thickness
boxWidthFrac = 0.50;   % box width as fraction of chord

m_wing_side = 11000;   % structural mass of one wing (kg)
m_engine    = 7700;    % installed engine + nacelle + pylon (kg)
yEngine     = 9.6;     % engine position from aircraft centreline (m)
m_fuel_side = 0;       % fuel in one wing (kg). 0 = no fuel relief
                       % (conservative up-bending case)

% --- Structure [assumptions] ---
eps_allow = 0.0045;    % allowable skin strain at ultimate load (0.45%)
E         = 70e9;      % effective Young's modulus of the CFRP wing box (Pa)
t_min     = 0.004;     % minimum equivalent skin thickness (m)

% --- Load cases ---
n_cases     = [1 2.5 3.75];   % 1g cruise, limit load, ultimate load
n_ult       = 3.75;
boeingTest  = 7.6;            % Boeing's measured tip deflection (m)

% --- Output ---
N           = 2000;           % spanwise points
saveFigures = true;
gifName     = 'wing_bending.gif';

%% =========================== 1. GEOMETRY =============================
g  = 9.81;
B  = b/2;                       % semi-span from centreline (m)
yf = fusW/2;                    % side of body (m)
Lam = deg2rad(sweep_deg);
yk = yf + yKink_fromSOB;

yc = linspace(0, B, N)';        % centreline-based coordinate (m)
chordFun = @(y) interp1([0 yf yk B], [c_root c_root c_kink c_tip], y);
tcFun    = @(y) interp1([0 yf yk B], [tc_root tc_root tc_kink tc_tip], y);
c  = chordFun(yc);

ex = yc >= yf;                  % exposed wing (outboard of fuselage)
y  = yc(ex);                    % exposed spanwise positions (m)
s  = (y - yf)/cos(Lam);         % distance along the swept beam axis (m)
cc = c(ex);
h  = boxDepthFrac * tcFun(y) .* cc;          % wing-box depth (m)
wbox = boxWidthFrac * cc * cos(Lam);         % box width normal to beam (m)

S_gross = 2*trapz(yc, c);
fprintf('Gross wing area   : %.0f m^2 (Wikipedia-derived: ~376 m^2)\n', S_gross);
fprintf('Aspect ratio      : %.2f (Wikipedia: 9.6)\n', b^2/S_gross);
fprintf('Beam length       : %.1f m along the swept wing\n', s(end));

% Everything the load function needs
P = struct('g',g,'W',MTOW,'B',B,'yc',yc,'c',c,'ex',ex,'y',y,'s',s, ...
           'cc',cc,'h',h,'Lam',Lam,'m_wing',m_wing_side, ...
           'm_fuel',m_fuel_side,'m_eng',m_engine,'sEng',(yEngine-yf)/cos(Lam));

%% ===================== 2-3. LOADS, SHEAR, MOMENT =====================
LU = wingLoads(n_ult, P);
fprintf('Root shear  at %.2fg : %.2f MN\n',   n_ult, LU.V(1)/1e6);
fprintf('Root moment at %.2fg : %.1f MN.m\n', n_ult, LU.M(1)/1e6);

%% ========================== 4. STIFFNESS =============================
% Method A: fully stressed wing -> curvature set by strain and depth only
kappaA_ult = 2*eps_allow ./ h;
wA_ult     = deflect(s, kappaA_ult);

% Method B: size skins for eps_allow at ultimate, with a minimum gauge
t_skin = max(LU.M ./ (E*eps_allow*wbox.*h), t_min);
EI     = E * wbox .* t_skin .* h.^2 / 2;
fprintf('Root equivalent skin thickness: %.1f mm\n', 1e3*t_skin(1));

%% ========================== 5. DEFLECTION ============================
nC = numel(n_cases);
wA = zeros(numel(s), nC);  wB = zeros(numel(s), nC);
for k = 1:nC
    Lk = wingLoads(n_cases(k), P);
    wB(:,k) = deflect(s, Lk.M ./ EI);
    wA(:,k) = (n_cases(k)/n_ult) * wA_ult;      % loads scale linearly with n
end

fprintf('\n=========== TIP DEFLECTION ===========\n');
fprintf(' Load case     Method A    Method B\n');
for k = 1:nC
    fprintf(' %4.2fg        %6.2f m    %6.2f m\n', n_cases(k), wA(end,k), wB(end,k));
end
fprintf(' Boeing test (3.75g): %.1f m\n', boeingTest);
fprintf(' Error vs test: A %+.0f%%, B %+.0f%%\n', ...
        100*(wA(end,end)/boeingTest-1), 100*(wB(end,end)/boeingTest-1));

%% ========================= 6. SENSITIVITY ============================
epsRange = linspace(0.0030, 0.0060, 25);
tipA_eps = zeros(size(epsRange));  tipB_eps = tipA_eps;
for k = 1:numel(epsRange)
    tipA_eps(k) = trapzTip(s, 2*epsRange(k)./h);
    ts  = max(LU.M ./ (E*epsRange(k)*wbox.*h), t_min);
    tipB_eps(k) = trapzTip(s, LU.M ./ (E*wbox.*ts.*h.^2/2));
end

thkScale = linspace(0.8, 1.2, 25);
tipB_thk = zeros(size(thkScale));
for k = 1:numel(thkScale)
    hk  = h*thkScale(k);
    ts  = max(LU.M ./ (E*eps_allow*wbox.*hk), t_min);
    tipB_thk(k) = trapzTip(s, LU.M ./ (E*wbox.*ts.*hk.^2/2));
end

%% =========================== 7. FIGURES ==============================
cols = [0.20 0.45 0.80; 0.95 0.55 0.10; 0.85 0.15 0.15];   % 1g, 2.5g, 3.75g

% ---------- Figure 1: planform and loads ----------
xqc = max(yc - yf, 0) * tan(Lam);            % quarter-chord line
xLE = xqc - 0.25*c;  xTE = xqc + 0.75*c;
f1 = figure('Color','w','Position',[50 50 1000 700]); lightTheme(f1);
subplot(2,1,1); hold on
fill([yc; flipud(yc)], [xLE; flipud(xTE)], [0.85 0.9 1], 'EdgeColor', [0.2 0.3 0.6]);
plot(yc, xqc, 'k--');  xline(yf, ':', 'Side of body');
plot(yEngine, interp1(yc, xLE, yEngine) - 2, 'k^', 'MarkerFaceColor', 'k');
set(gca, 'YDir', 'reverse'); axis equal; grid on; box on
xlabel('Spanwise position from centreline (m)'); ylabel('Chordwise (m)');
title(sprintf('787-8 half-wing planform model (S = %.0f m^2, sweep %d°)', S_gross, sweep_deg));
subplot(2,1,2); hold on
plot(y, LU.lift/1e3, 'LineWidth', 2);
plot(y, -LU.wingWt/1e3, 'LineWidth', 2);
plot(y, LU.lift/1e3 - LU.wingWt/1e3 - LU.fuelWt/1e3, 'k', 'LineWidth', 2);
xlabel('Spanwise position from centreline (m)'); ylabel('Load (kN/m)');
legend({'Lift (Schrenk)','Wing weight','Net load'}, 'Location', 'northeast');
title(sprintf('Spanwise loads at ultimate load (n = %.2f)', n_ult));
grid on; box on

% ---------- Figure 2: shear and bending moment ----------
f2 = figure('Color','w','Position',[80 80 1000 600]); lightTheme(f2);
subplot(2,1,1);
plot(s, LU.V/1e6, 'LineWidth', 2); grid on; box on
ylabel('Shear force (MN)'); title(sprintf('Shear force and bending moment at %.2fg', n_ult));
subplot(2,1,2);
plot(s, LU.M/1e6, 'LineWidth', 2); grid on; box on
xlabel('Distance along wing from side of body (m)'); ylabel('Bending moment (MN·m)');

% ---------- Figure 3: deflection curves ----------
f3 = figure('Color','w','Position',[110 110 950 520]); lightTheme(f3); hold on
hB = gobjects(nC,1);
for k = 1:nC
    hB(k) = plot(y, wB(:,k), '-', 'Color', cols(k,:), 'LineWidth', 2.5);
end
hA = plot(y, wA(:,end), '--', 'Color', cols(end,:), 'LineWidth', 1.5);
hT = plot(y(end), boeingTest, 'kp', 'MarkerSize', 16, 'MarkerFaceColor', 'y');
text(y(end)-1, boeingTest-0.2, 'Boeing test: 7.6 m', 'HorizontalAlignment', 'right', 'FontSize', 12);
ylim([0 boeingTest+0.8]);
xlabel('Spanwise position from centreline (m)'); ylabel('Upward deflection (m)');
title('787-8 wing deflection: model vs Boeing''s 2010 ultimate-load test');
legend([hB; hA; hT], {'1g cruise','2.5g limit load','3.75g ultimate load', ...
       '3.75g (Method A: fully stressed)','Boeing measured'}, 'Location', 'northwest');
grid on; box on; set(gca, 'FontSize', 12);

% ---------- Figure 4: 3D bent wing ----------
f4 = figure('Color','w','Position',[140 140 1100 600]); lightTheme(f4); hold on
yE = y;  xLEe = xLE(ex);  xTEe = xTE(ex);
for k = 1:nC
    Z = [wB(:,k) wB(:,k)];
    surf([xLEe xTEe], [ yE  yE], Z, 'FaceColor', cols(k,:), 'EdgeColor', 'none', 'FaceAlpha', 0.75);
    surf([xLEe xTEe], [-yE -yE], Z, 'FaceColor', cols(k,:), 'EdgeColor', 'none', 'FaceAlpha', 0.75);
end
Z0 = zeros(numel(yE), 2);
mesh([xLEe xTEe], [ yE  yE], Z0, 'EdgeColor', [0.5 0.5 0.5], 'FaceColor', 'none');
mesh([xLEe xTEe], [-yE -yE], Z0, 'EdgeColor', [0.5 0.5 0.5], 'FaceColor', 'none');
drawFuselage(yf);
axis equal; grid on; view(-65, 10);                % looking from the front-left
ylim([-B-1 B+1]); zlim([-3.5 boeingTest+1]);
xlabel('x (m)'); ylabel('Span (m)'); zlabel('Deflection (m)');
title('787-8 wing shape: unloaded (grey), 1g (blue), 2.5g (orange), 3.75g (red)');
camlight headlight; lighting gouraud; set(gca, 'FontSize', 11);

% ---------- Figure 5: sensitivity ----------
f5 = figure('Color','w','Position',[170 170 1000 450]); lightTheme(f5);
subplot(1,2,1); hold on
plot(100*epsRange, tipA_eps, '--', 'LineWidth', 2);
plot(100*epsRange, tipB_eps, '-', 'LineWidth', 2);
yline(boeingTest, 'k:', 'Boeing 7.6 m', 'LineWidth', 1.5);
xlabel('Allowable skin strain at ultimate (%)'); ylabel('Tip deflection at 3.75g (m)');
legend({'Method A','Method B'}, 'Location', 'northwest');
title('Effect of strain allowable'); grid on; box on
subplot(1,2,2); hold on
plot(100*thkScale, tipB_thk, '-', 'LineWidth', 2);
yline(boeingTest, 'k:', 'Boeing 7.6 m', 'LineWidth', 1.5);
xlabel('Wing thickness (% of baseline)'); ylabel('Tip deflection at 3.75g (m)');
title('Thinner wing = more bending'); grid on; box on

if saveFigures
    saveFig(f1, 'fig1_planform_loads.png');
    saveFig(f2, 'fig2_shear_moment.png');
    saveFig(f3, 'fig3_deflection.png');
    saveFig(f4, 'fig4_wing_3d.png');
    saveFig(f5, 'fig5_sensitivity.png');
end

% ---------- GIF: wing bending from 0 to 3.75g (front view) ----------
f6 = figure('Color','w','Position',[200 200 1000 420]); lightTheme(f6);
if exist(gifName, 'file'), delete(gifName); end
nAnim = [linspace(0, n_ult, 40), n_ult*ones(1,10)];
for k = 1:numel(nAnim)
    clf(f6);
    ax6 = axes(f6); hold(ax6, 'on');          % draw only into this figure
    wk = wB(:,end) * nAnim(k)/n_ult;
    plot(ax6, [-flipud(y); y], [flipud(wk); wk], 'Color', cols(end,:), 'LineWidth', 4);
    plot(ax6, [-yf yf], [0 0], 'k', 'LineWidth', 8);
    plot(ax6, [-B B], [0 0], ':', 'Color', [0.6 0.6 0.6]);
    axis(ax6, 'equal'); xlim(ax6, [-B-1 B+1]); ylim(ax6, [-1 boeingTest+1.5]);
    grid(ax6, 'on'); box(ax6, 'on');
    title(ax6, sprintf('Load factor n = %.2fg   |   tip deflection = %.2f m', nAnim(k), wk(end)), 'FontSize', 13);
    xlabel(ax6, 'Span (m)'); ylabel(ax6, 'Deflection (m)');
    ax6.Toolbar = [];
    drawnow;
    exportgraphics(f6, gifName, 'Append', true);
end
fprintf('\nSaved figures and %s\n', gifName);

%% =====================================================================
%%                          LOCAL FUNCTIONS
%% =====================================================================
function out = wingLoads(n, P)
% Net loads, shear and bending moment along the swept beam at load factor n.
    Lh   = n*P.W*P.g/2;                                     % lift on one side (N)
    le   = 4*Lh/(pi*P.B) * sqrt(max(1 - (P.yc/P.B).^2, 0)); % elliptical
    lp   = Lh * P.c / trapz(P.yc, P.c);                     % planform-shaped
    lift = 0.5*(le + lp);                                   % Schrenk (N/m)
    lift = lift(P.ex);
    wingWt = n*P.g*P.m_wing * P.cc / trapz(P.y, P.cc);
    vol    = P.cc .* P.h;
    fuelWt = n*P.g*P.m_fuel * vol / trapz(P.y, vol);
    q = (lift - wingWt - fuelWt) * cos(P.Lam);              % per metre of beam
    V = -flipud(cumtrapz(flipud(P.s), flipud(q)));          % shear (tip -> root)
    V = V - n*P.g*P.m_eng*(P.s < P.sEng);                   % engine point load
    M = -flipud(cumtrapz(flipud(P.s), flipud(V)));          % bending moment
    out = struct('lift',lift,'wingWt',wingWt,'fuelWt',fuelWt,'V',V,'M',M);
end

function w = deflect(s, kappa)
% Integrate curvature twice from the clamped root (w = 0, slope = 0).
    w = cumtrapz(s, cumtrapz(s, kappa));
end

function tip = trapzTip(s, kappa)
    w = deflect(s, kappa);  tip = w(end);
end

function drawFuselage(r)
    [Yc, Zc, Xc] = cylinder(r, 40);
    surf(Xc*56.7 - 22, Yc, Zc - 1.0, 'FaceColor', [0.8 0.8 0.8], ...
         'EdgeColor', 'none', 'FaceAlpha', 0.35);
end

function lightTheme(f)
    try, theme(f, 'light'); catch, end
end

function saveFig(f, name)
    if ~isgraphics(f)
        warning('%s skipped: the figure window was closed before saving.', name);
        return
    end
    set(findall(f, 'Type', 'axes'), 'Toolbar', []);
    exportgraphics(f, name, 'Resolution', 200);
end
