set t
/1*12/
    it
/i1*i12/
    branch  (it,t) 'time periods in each iteration'
/i1.(1*4),
i2.(2*5),
i3.(3*6),
i4.(4*7),
i5.(5*8),
i6.(6*9),
i7.(7*10),
i8.(8*11),
i9.(9*12),
i10.(10*12),
i11.(11*12),
i12.(12)/
    opt_time(t,it)    'optimization time'
    tfirst(t)
    first_t_it(t,it)  'first hour in each iteration'
    Second_t_it(t,it) 'second hour in each iteration'
   
;

tfirst(t) = yes$(ord(t) eq 1);


set
        p        power plants           /
P1_NUC
P2_OCGT
P3_LIG
P4_CCGT
P5_Wind
/
        c(p)            conventional plants
        r(p)            RES plants
;
alias(it,iitt)
alias(t,tt)
;
PARAMETERS
        genup                   all information about generation
        timeup                  all information about time-series
        cap(p)                  capacities of PPs
        vc(p)                   variable production costs of PPs
        demand_forecast(t,it)      demand forecasted
        demand_real (t,it)         actual demand
        af_wind(t)              availability of wind
        SU_costs(p)             start_up_costs
;

$onecho > input_data.tmp
        par=genup            rng=Power!A1                  cdim=1 rdim=1
        par=timeup           rng=Time!A1                   cdim=1 rdim=1
 
$offecho

$onUNDF
$call   gdxxrw I=try_data.xlsx O=try_data.gdx cmerge=1 @input_data.tmp
$gdxin  try_data.gdx
$load   genup, timeup
$offUNDF
 
*execute_unload "check.gdx";
*$stop

c(p)    = genup(p,'type')= 1 ;
r(p)    = genup(p,'type')= 2 ;
cap(p)  = genup(p,'capacity');
vc(p)                   = genup(p,'vc')      ;
demand_forecast(t,it)   = timeup(t,'Demand_forecast')  ;
demand_real(t,it)       = timeup(t,'actual_demand')  ;
af_wind(t)              = timeup(t,'Wind_af') ;
SU_costs(p)             = genup(p,'Start_up_costs');

*execute_unload "check.gdx";
*$stop
   
Parameters
count

price(t,it)

cap_stor              
power_turb               
power_pump               
loss_stor
g_min            
;


count = 1;

cap_stor    = 200;   
power_turb  = 100;     
power_pump  = 50;    
loss_stor   = 0.2;
g_min = 0.3;



variable
        TC
;

positive variable
G(p,t,it)
report_gen(p,t,it,*)
report_stor(t,it,*)

level_stor(t,it)         
charge_stor(t,it)       
gen_stor(t,it)           




P_ON(p,t,it) running capacity    [MW]
SU(p,t,it)   start-up activity   [MW]
;

equations

OBJECTIVE
ENERGY_BALANCE
MAX_GENERATION
MAX_GENERATION_WIND

stor_level_def
Stor_level_RH
stor_level_cap        
stor_gen_cap_level      
stor_charge_cap

*stor_level_iteration


min_generation
max_online
startup
;

OBJECTIVE..                                          TC =e= sum((p,t,iitt)$opt_time(t,iitt), G(p,t,iitt)*vc(p))
                                                            + sum((c,t,iitt), Su(c,t,iitt) * SU_costs(c));

ENERGY_BALANCE(t,iitt)$opt_time(t,iitt)..            demand_forecast(t,iitt) =e= sum(p, G(p,t,iitt))
                                                                                + gen_stor(t,iitt) - charge_stor(t,iitt);
    

MAX_GENERATION(c,t,iitt)..                     G(c,t,iitt)        =l= P_On(c,t,iitt);
*                                                                    cap(c)
MAX_GENERATION_WIND(r,t,iitt)..                G(r,t,iitt)        =l= cap(r)*af_wind(t);


*******************Storage*******************
stor_level_def(t,iitt)$opt_time(t,iitt)..            level_stor(t+1,iitt)  =e= level_stor(t,iitt) + charge_stor(t+1,iitt)*(1-loss_stor) - gen_stor(t+1,iitt);
Stor_level_RH(t,iitt)$opt_time(t,iitt)..             level_stor(t,iitt)$tfirst(t)    =e= 0;
stor_level_cap(t,iitt)..            level_stor(t,iitt)  =l= cap_stor  ;
stor_gen_cap_level(t,iitt)..        gen_stor(t,iitt)    =l= power_turb;
stor_charge_cap(t,iitt)..           charge_stor(t,iitt) =l= power_pump;


*stor_level_iteration(t,iitt,first_t_it,second_t_it)..   level_stor(t,iitt)$first_t_it  =e= level_stor(t,iitt-1)$second_t_it;

******************Start_Up_Activity**********

min_generation(c,t,iitt)..                             P_on(c,t,iitt) * g_min                   =l= G(c,t,iitt)           ;
max_online(c,t,iitt)..                                 G(c,t,iitt)                              =l= cap(c)        ;
startup(c,t,iitt)..                                    P_on(c,t,iitt) - P_on(c,t-1,iitt)        =l= SU(c,t,iitt)          ;


model try /all/;

*execute_unload "check.gdx";
*$stop



Loop (it,

opt_time(t,it)  =    branch (it,t);
 
if (ord(it) eq count,
demand_forecast(t,it)$(ord(t)=count)= demand_real(t,it);
demand_forecast(t,it)$(ord(t) le count -1) = 0;

else
demand_forecast(t,it)= demand_forecast(t,it);
);

first_t_it(t,it) = count$(ord(t)=count);
Second_t_it(t,it)= count$(ord(t)=count+1);


solve try using lp minimizing TC ;

report_gen.l(p,t,it,'generation') = G.l(p,t,it);

*G.fx(p,t,it)$(ord(t)= ord(it) and opt_time(t,it)) = G.l(p,t,it);

report_stor.l(t,it,'storage') = level_stor.l(t,it);
level_stor.fx(t,it) = level_stor.l(t,it);


P_on.fx(c,t,it) = P_on.l(c,t,it);

Price(t,it)$opt_time(t,it) = - Energy_balance.m(t,it);

count = count+1;


);




execute_unload "check.gdx";
$stop
