Sets p power plants
 /

P1_NUC
P2_OCGT
P3_LIG
P4_CCGT
P5_Wind
P6_CCGT
P7_OCGT
P8_Wind
P9_NUC
P10_OCGT
P11_LIG
P12_CCGT
P13_Wind
P14_CCGT
P15_OCGT
P16_Wind
/
 
    t time periods
/1*48/
    
    it iterations
/iteration1*iteration12/
*until iteration48

  Wind_it  wind availability forcasts
/Wind_availability1 * Wind_availability12/
* until Wind_availability48


*execute_unload "check.gdx";
*$stop


;
SETS
    c(p)            conventional plants
    r(p)            RES plants

;

PARAMETERS
        Generation               all information about generation
        cap(p)                   capacities of PPs
        vc(p)                    variable production costs of PPs
        demand(t,it)             demand in in every periode      
        af_wind(t,wind_it)       availability factor for wind feed-in in every periode
        
        PH(t)                    prediction horizon
        CH(t)                    controll horizon
        
        SU_costs(p,t,it)
        SH_costs(p,t,it)
;


$onecho > input_data.tmp
         par=Generation    rng=Generation!A1      cdim=1 rdim=1
         par=demand        rng=Timing_RH!A1       cdim=1 rdim=1
         par=af_wind       Rng=Timing_wind!A1     cdim=1 rdim=1
$offecho

$onUNDF
$call   gdxxrw I=Input_RH.xlsx O=Input_RH.gdx cmerge=1 @input_data.tmp
$gdxin  Input_RH.gdx
$load   Generation, demand, af_wind
$offUNDF

*execute_unload "check.gdx";
*$stop

c(p)    = Generation(p,'Type')= 1;
r(p)    = Generation(p,'Type')= 2 ;
cap(p)  = Generation(p,'Capacity');
vc(p)   = Generation(p,'VC');
SU_costs(p,t,it) = Generation(p,'start up costs');
SH_costs(p,t,it)= Generation(p,'shut down costs');


*execute_unload "check.gdx";
*$stop

variable TC
;
positive variable
G(p,t,it);

binary variable
Start_up(p,t,it)
Shut_down(p,t,it);

*execute_unload "check.gdx";
*$stop

equations

OBJECTIVE
ENERGY_BALANCE
MAX_GENERATION
MAX_GENERATION_WIND

;
Objective..  TC =e= sum((p,t,it),G(p,t,it) * vc(p))
                +   sum((p,t,it),Start_up(p,t,it)* SU_costs(p,t,it))
                +   sum((p,t,it),SH_costs(p,t,it)*SH_costs(p,t,it) )
;

Energy_balance(t,it)..  demand(t,it) =e= sum((p), G(p,t,it));

Max_Generation(c,t,it)..  G(c,t,it) =l= cap(c)* Start_up(c,t,it);

Max_Generation_wind(r,t,it,wind_it)..    G(r,t,it) =l= cap(r)* af_wind(t,wind_it);

model RH /all/;
*af_wind(t,wind_it)= 0;



solve RH using MIP minimizing TC;
*loop((t),
*tph(t) = tph(t+1);
*)

*execute_unload "check.gdx";
*$stop


Parameters

price(t,it);

Price(t,it) = - Energy_balance.m(t,it);


execute_unload "check.gdx";
$stop












