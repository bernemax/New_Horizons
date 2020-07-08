set 
it
/1*10/
t
/1*12/

i(it,t)

opt_time(it,t)    'optimization time'
;
i(it,t)$(ord(t) ge ord(it) and ord(t) le (ord(it)+3) )= yes;
     
set
first_t_it(t,it)  'first hour in each iteration'
Second_t_it(t,it) 
*execute_unload "check.gdx";
*$stop
;
set
        p        power plants
/
P1_NUC
P2_OCGT
P3_LIG
P4_CCGT
P5_Wind
/
        c(p)            conventional plants
        r(p)            RES plants
;

alias (t,tt);
alias (it,iitt);

PARAMETERS
        genup                   all information about generation
        timeup                  all information about time-series
        cap(p)                  capacities of PPs
        vc(p)                   variable production costs of PPs
        demand_forecast(t,it)      demand forecasted
        demand_real (t,it)         actual demand
        af_wind(t)              availability of wind
        SU_costs(p)             start_up_costs
        count
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
count                   = 1;



*execute_unload "check.gdx";
*$stop

loop(it,

  if (ord(it) eq count,
        demand_forecast(t,it)$(ord(t)=count)= demand_real(t,it);
        demand_forecast(t,it)$(ord(t) le count -1) = demand_real(t,it);
        demand_forecast(t,it)$(ord(t) gt count +3) = 0;
        demand_real(t,it)$(ord(t) gt count )       = 0;
        demand_forecast(t,it)$(ord(t) < ord(it)) = 0;
    else
        demand_forecast(t,it)= demand_forecast(t,it);
        );
    count= count +1;
);
*execute_unload "check.gdx";
*$stop   
Parameters
low
high
price(t,it)
g_min
count
online(p,t,it)

test;

g_min = 0.3;
low =1;
High =4;
count=1;
test = 1;
online(p,t,it)=0;

variable
        TC
;

positive variable
G(p,t,it)

P_ON(p,t,it) running capacity    [MW]
SU(p,t,it)   start-up activity   [MW]
;

equations

OBJECTIVE
ENERGY_BALANCE
MAX_GENERATION
MAX_GENERATION_WIND

min_generation
max_online
startup
*start_up_rest
;

OBJECTIVE..                                          TC =e= sum((p,t,iitt)$opt_time(iitt,t) , G(p,t,iitt)*vc(p))
                                                          + sum((c,t,iitt), Su(c,t,iitt) * SU_costs(c));

ENERGY_BALANCE(t,iitt)$opt_time(iitt,t)..                            demand_forecast(t,iitt) =e= sum(p, G(p,t,iitt));
*                                                                               + gen_stor(t,iitt) - charge_stor(t,iitt);
    

MAX_GENERATION(c,t,iitt)..                     G(c,t,iitt)        =l= P_On(c,t,iitt);
*                                                                    cap(c)
MAX_GENERATION_WIND(r,t,iitt)..                G(r,t,iitt)        =l= cap(r)*af_wind(t);

******************Start_Up_Activity**********

min_generation(c,t,iitt)..                             P_on(c,t,iitt) * g_min                   =l= G(c,t,iitt)           ;
max_online(c,t,iitt)..                                 G(c,t,iitt)                             =l= cap(c)                ;
startup(c,t,iitt)..                                    P_on(c,t,iitt) - P_on(c,t-1,iitt)$(ord(t) ge ord (iitt))   =l= SU(c,t,iitt)          ;
            
model try /all/;

*execute_unload "check.gdx";
*$stop

$ontext
Loop (branch,

*[time definitions - which hours to solve]
1) if branch is 1: hours 1-3 in model
2) now(t) t
3) tomorrow(t) t+1
4) d+1(t) t+1
*solve model
*export data
*);
$offtext


loop (it,

opt_time(it,t) = i(it,t) ;

first_t_it(t,it) = count$(ord(t)=count);
Second_t_it(t,it)= count$(ord(t)=count+1);

*$ontext    
    if (ord(it) > 1,
        p_on.l(c,t,it)$(first_t_it(t,it)) = online(c,t,it-1)$(ord(t) = (ord(it)));
        solve try using lp minimizing tc;
        online(c,t,it)= p_on.l(c,t,it);
    
    else
    
        solve try using lp minimizing tc;
        online(c,t,it)= p_on.l(c,t,it) ;  
    );
*$offtext    
    
    price(t,it)= -energy_balance.m(t,it);
                   
count= count +1;      
    
);
execute_unload "check.gdx";
$stop   
    
    




*put_utility 'gdxout' / 'output_v20_EVP\v20_STOCH_FINV_Scen_' UncPar.tl:0;
*execute_unload  ;
