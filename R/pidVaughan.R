library(ggplot2)
library(FLCore)
library(reshape)

if (FALSE){
#yq = []; yv = []; ya = []; yj = [];

q = .3; v = .1; a = .2; 

setPt = 0.0
amax  = 0.3 
h     = 1.0

dt = .01; r = seq(0,12,dt);
lim = 12;

y=array(0,dim     =c(4,length(r)),
          dimnames=list(quant=c("q","v","a","j"),
                        t    =round(r/dt) + 1))

for (i in ac(round(r/dt) + 1)){
  if (q < setPt){
    q = setPt - q
    v = -v
    a = -a
    }
  
  if (v > 0)
    F = h
  else
    F = setPt
  
  if (abs(v*dt) < abs(F - q)){
    p = v/(F - q);
    if (v > 0)
      j = -(amax + a)*p
    else 
      j = -(3*a + 2*v*p)*p
  }
  
  y[4,i] = j;
  a = a + j*dt; y[3, i] = a;
  v = v + a*dt; y[2, i] = v;
  q = q + v*dt; y[1, i] = q;
  }

nms=c('position', 'velocity', 'acceleration', 'jerk')
names(nms)=substr(nms,1,1)
  
dat=melt(cbind(r=r,data.frame(t(y))),id="r")
dat$variable=nms[dat$variable]

ggplot(dat)+
  geom_line(aes(r,value,colour=variable),size=2)+
  theme_bw()+
  theme(legend.position="bottom")

# % figure(1); clf; grid on; box on; axis auto; hold on
# co = get(gca, 'ColorOrder'); ls = {'-' '--' ':' '-.'};
# plot(r, y(1, , r, y(2, , r, y(3, , r, y(4, , 'linewidth', 2);
#     % for j = 1:4
#     %     ls{j}
#     %     line(range, y(j, , 'linewidth', 2, 'Color', co(j, , 'LineStyle', '--')
#     % end
#     legend('position', 'velocity', 'acceleration', 'jerk'); grid on
#     text(1, .8, 'Sinusoid out', 'Color', 'b', 'FontSize', 18)
#     text(5, .67, 'Gaussian return', 'Color', 'b', 'FontSize', 18)
#     % text(bangt, bangy + h/20, '| Bang at apogee', 'Color', 'b', 'FontSize', 18)
#     % text(dampt, dampy + h/10, '| Damping onset at j = 0', 'Color', [0 .75 .75], 'FontSize', 18)
#     set(title('Maneuver: avoid limit at L = 1m and return to setpoint at S = 0m'), 'fontsize', 20)
#     set(gca, 'xgrid', 'on', 'ygrid', 'on', 'ytick', -20:.5:20, 'xlim', [0 lim], 'ylim', [-.5 1])
#     % print -dwin SinGau.eps
}