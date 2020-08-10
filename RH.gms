set t scheduling horizon /1*216/

    it number of iterations /1*9/
    
    day(t,it)
    
    CH(t,it) calculation horizon
    
    PH(t,it) prediction horizon
    
    FH(t,it) fixed horizon
    
   
;

*    CH(t,it)$(ord(t) ge ord(it) + (23*(ord(it)-1))$(ord(it)>1) and ord(t) le (ord(it)*24) ) = yes;
*    CH(t,it)$(ord(t) eq  (24*(ord(it)-1))$(ord(it)>1)) = yes;
    

*    PH(t,it)$(ord(t) ge ord(it) + (23*(ord(it)-1))$(ord(it)>1) and ord(t) le (ord(it)*24) ) = yes;
    
    
*    FH(t,it)$(ord(t) eq ord(it) * 24 ) = yes;
*    FH(t,it)$(ord(t) eq  (24*(ord(it)-1))$(ord(it)>1)) = yes;
;

    
set
    p power plants
/
P1_NUC
P2_OCGT
P3_LIG
P4_CCGT
P5_Wind
/
    c(p) conventional plants
    
    r(p) RES plants
;

alias (t,tt);
alias (it,iitt);

PARAMETERS
        genup                   all information about generation
        timeup                  all information about time-series
        cap(p)                  capacities of PPs
        vc(p)                   variable production costs of PPs
        demand_forecast(t,it)   demand forecasted
        af_wind_forecast(t,it)  forecasted availability of wind
        af_wind_actual(t,it)    actual availability of wind
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
 

c(p)                    = genup(p,'type')= 1 ;
r(p)                    = genup(p,'type')= 2 ;
cap(p)                  = genup(p,'capacity');
vc(p)                   = genup(p,'vc')      ;
demand_forecast(t,it)   = timeup(t,'Demand_forecast')  ;
af_wind_forecast(t,it)  = timeup(t,'Wind_af_forecast') ;
af_wind_actual(t,it)    = timeup(t,'wind_af_actual');
SU_costs(p)             = genup(p,'Start_up_costs');
count                   = 1;

$ontext
loop(it,

    if (ord(it),
        demand_forecast(t,it)$(CH(t,it))= demand_forecast(t,it);
    
     else
        demand_forecast(t,it) = 0;
        );
        
  
);
$Offtext


Parameters
price(t,it)
g_min
count
online(p,t,it)
generation(p,t,it)

test;

g_min = 0.3;
count = 1;
test  = 1;
*online(p,t,it)=0;

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
ENERGY_BALANCE(t,it)
MAX_GENERATION(p,t,it)
MAX_GENERATION_WIND(p,t,it)

min_generation(p,t,it)
max_online(p,t,it)
startup(p,t,it)

;

OBJECTIVE..                                                 TC =e= sum((p,CH), G(p,Ch)*vc(p))
                                                             + sum((c,Ch), Su(c,CH) * SU_costs(c));

ENERGY_BALANCE(CH)..                                    demand_forecast(CH) =e= sum(p, G(p,CH));
*                                                                               + gen_stor(t,iitt) - charge_stor(t,iitt);
    

MAX_GENERATION(c,CH)..                                  G(c,CH)        =l= P_On(c,CH);
*                                                                    cap(c)
MAX_GENERATION_WIND(r,CH)..                             G(r,CH)        =l= cap(r)*af_wind_forecast(CH);

******************Start_Up_Activity**********

min_generation(c,CH)..                                  P_on(c,CH) * g_min                        =l= G(c,CH)               ;
max_online(c,CH)..                                      G(c,CH)                                   =l= cap(c)                ;
startup(c,t,iitt)$(PH(t,iitt))..                        P_on(c,t,iitt) - P_on(c,t-1,iitt)         =l= SU(c,t,iitt)          ;
            
model try /all/;

loop(it ,

**************************************time definition in loop**************************************
    CH(t,it)$(ord(t) ge ord(it) + (23*(ord(it)-1))$(ord(it)>1) and ord(t) le (ord(it)*24) ) = yes;
    CH(t,it)$(ord(t) eq  (24*(ord(it)-1))$(ord(it)>1)) = yes;
    

    PH(t,it)$(ord(t) ge ord(it) + (23*(ord(it)-1))$(ord(it)>1) and ord(t) le (ord(it)*24) ) = yes;
    
    
    FH(t,it)$(ord(t) eq ord(it) * 24 ) = yes;
    FH(t,it)$(ord(t) eq  (24*(ord(it)-1))$(ord(it)>1)) = yes;

**************************************demand definition over time************************************
    if (ord(it) eq count,
        demand_forecast(t,it)$(CH(t,it))= demand_forecast(t,it);
    
     else
        demand_forecast(t,it) = 0;
    ); 

solve try using lp minimizing tc;
  
P_on.fx(p,t,it+1)$(FH(t,it))= P_on.l(p,t,it)$(FH(t,it));

price(t,it)= -energy_balance.m(t,it);

count=count+1;
);
execute_unload "check.gdx";  
$stop







****************************************Results of the last skype call**********************************
$ontext
*[time definitions - which hours to solve]
1) if branch is 1: hours 1-3 in model
2) now(t) t
3) tomorrow(t) t+1
4) d+1(t) t+1
*solve model
*export data
*);
$offtext


*put_utility 'gdxout' / 'output_v20_EVP\v20_STOCH_FINV_Scen_' UncPar.tl:0;
*execute_unload  ;
