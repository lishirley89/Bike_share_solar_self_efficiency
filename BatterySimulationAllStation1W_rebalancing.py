#simulate with calculated demand from calcDemand.py and pv from calcPV.py
#consider rebalancing, replace battery during rebalancing and battery level lower than 40% (RebalanceRefill)
#battery stop discharging at 20% (BatteryCutout)
from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir("/users/shirley/Documents/Research/Chicago/Battery/MyBatteryModel")

times = pd.date_range(start='2019-01-01', end='2019-12-31 23:00:00', freq='H')
#load demand if needed
demand_all = pd.read_csv('results/demand2019.csv',index_col=0)
demand_all =demand_all.set_index(times)

#load pv
pv_chicago = pd.read_csv('data/pv_chicago_year.csv',index_col=0)
pv_chicago =pv_chicago.set_index(times)

#load rebalancing record for all stations
rebalance_all = pd.read_csv('results/reb_hour2019.csv',index_col=0)
rebalance_all=rebalance_all.set_index(times)

# notice battery stop discharging at 20%
param_tech = {'BatteryCapacity': 1.2,
              'BatteryEfficiency': 0.8, #attention!
              'InverterEfficiency': 0.96,
              'timestep': 1,
              'MaxPower': 0.114,
              'BatteryCutout': 0.2,
              'RebalanceRefill': 0.4
             }



#########################################################
#         dispatch function, off-grid, start level full capacity
#########################################################

def dispatch_max_sc(pv, demand, rebalance,param, return_series=False):
    """ Self consumption maximization pv + battery dispatch algorithm.
    The dispatch of the storage capacity is performed in such a way to maximize self-consumption:
    the battery is charged when the PV power is higher than the load and as long as it is not fully charged.
    It is discharged as soon as the PV power is lower than the load and as long as it is not fully discharged.

    Arguments:
        pv (pd.Series): Vector of PV generation, in kW DC (i.e. before the inverter)
        demand (pd.Series): Vector of household consumption, kW
        param (dict): Dictionary with the simulation parameters:
                timestep (float): Simulation time step (in hours)
                BatteryCapacity: Available battery capacity (i.e. only the the available DOD), kWh
                BatteryEfficiency: Battery round-trip efficiency, -
                InverterEfficiency: Inverter efficiency, -
                MaxPower: Maximum battery charging or discharging powers (assumed to be equal), kW
                BatteryCutout: cut out point for the power inverter (convert DC to AC), battery stop discharging
        return_series(bool): if True then the return will be a dictionary of series. Otherwise it will be a dictionary of ndarrays.
                        It is reccommended to return ndarrays if speed is an issue (e.g. for batch runs).
    Returns:
        dict: Dictionary of Time series

    """

    bat_size_e_adj = param['BatteryCapacity']
    bat_size_p_adj = param['MaxPower']
    n_bat = param['BatteryEfficiency']
    n_inv = param['InverterEfficiency']
    timestep = param['timestep']
    bat_cut_out = param['BatteryCutout']
    bat_reb_refill = param['RebalanceRefill']
    # We work with np.ndarrays as they are much faster than pd.Series
    Nsteps = len(pv)
    LevelOfCharge = np.zeros(Nsteps)
    pv2store = np.zeros(Nsteps)
    store2inv = np.zeros(Nsteps)
    # record if the battery is replaced with full
    bat_rep = np.zeros(Nsteps)
    # record battery leftover before replacement
    #bat_leftover = np.zeros(Nsteps)
    # record battery level below cut out point
    bat_cutout=np.zeros(Nsteps)

    #Load served by PV
    pv2inv = np.minimum(pv, demand / n_inv)  # DC direct self-consumption

    #Residual load
    res_load = (demand - pv2inv * n_inv)  # AC
    inv2load = pv2inv * n_inv  # AC

    #Excess PV
    res_pv = np.maximum(pv - demand/n_inv, 0)  # DC

    #PV to storage after eff losses
    pv2inv = pv2inv.values

    #########################################################
    # charge level start with full capacity
    #########################################################
    #first timestep = 0
    LevelOfCharge[0] = bat_size_e_adj   # DC
    for i in range(1,Nsteps):
        #PV to storage
        if LevelOfCharge[i-1] >= bat_size_e_adj:  # if battery is full
                pv2store[i] = 0
        else: #if battery is not full
            if LevelOfCharge[i-1] + res_pv[i] * n_bat * timestep > bat_size_e_adj:  # if battery will be full after putting excess
                pv2store[i] = min((bat_size_e_adj - LevelOfCharge[i-1]) / timestep, bat_size_p_adj)
            else:
                pv2store[i] = min(res_pv[i] * n_bat, bat_size_p_adj)

        # Storage to load
        store2inv[i] = min(bat_size_p_adj,  # DC
                           res_load[i] / n_inv,
                           LevelOfCharge[i-1] / timestep)

        #SOC
        LevelOfCharge[i] = min(LevelOfCharge[i - 1] - (store2inv[i] - pv2store[i]) * timestep,  # DC
                               bat_size_e_adj)
        # if battery is below rebalance refill point after putting, check if there is rebalance
        if LevelOfCharge[i - 1] + res_pv[i] * n_bat * timestep < bat_size_e_adj * bat_reb_refill:
            if rebalance[i] == 1:
                LevelOfCharge[i] = bat_size_e_adj  # assume a full battery will replace during rebalance
                bat_rep[i] = 1  # 1 means battery replaced
                # bat_leftover[i - 1] = max(LevelOfCharge[i - 1] + res_pv[i] * n_bat * timestep, 0)  # battery leftover

            else:
                LevelOfCharge[i] = min(LevelOfCharge[i - 1] - (store2inv[i] - pv2store[i]) * timestep,  # DC
                                       bat_size_e_adj)

                # if battery is below cut out point after putting excess, station down, no consumption, save energy in battery
                if LevelOfCharge[i - 1] + res_pv[i] * n_bat * timestep < bat_size_e_adj * bat_cut_out:
                    pv2inv[i] = 0  # save all pv in battery
                    pv2store[i] = min(pv[i] * n_bat, bat_size_p_adj)
                    store2inv[i] = 0  # stop discharging
                    inv2load[i] = 0  # station down, stop powering
                    res_pv[i] = 0  # no pv left
                    bat_cutout[i] = 1  # 1 means battery is shut down
                    LevelOfCharge[i] = min(LevelOfCharge[i - 1] + pv2store[i] * timestep,  # DC
                                           bat_size_e_adj)

    pv2store=pv2store/n_bat
    #pv2inv = pv2inv + res_pv - pv2store
    ### pv waste
    pv2waste = np.maximum(pv - pv2inv - pv2store,0)


    inv2load = inv2load + store2inv * n_inv  # AC

    out = {'pv2inv': pv2inv,
            'res_pv': res_pv,
            'pv2store': pv2store,
            'inv2load': inv2load,
            'store2inv': store2inv,
            'LevelOfCharge': LevelOfCharge,
            'BatteryReplace': bat_rep,
            #'BatteryLeftover': bat_leftover,
            'BatteryCutout': bat_cutout,
            'pv2waste': pv2waste
            }
    if not return_series:
        out_pd = {}
        for k, v in out.items():  # Create dictionary of pandas series with same index as the input pv
            out_pd[k] = pd.Series(v, index=pv.index)
        out = out_pd
    return out

#### battery replacement #####
BatteryReplaceAll_df=pd.DataFrame() # to hold all batteryreplace values for each station
#ConsumptionAll_df=pd.DataFrame() # to hold all Consumption/inv2load values for each station
#BatteryLeftoverAll_df=pd.DataFrame() # to hold all batteryleftover values for each station
BatteryCutoutAll_df=pd.DataFrame() # to hold all batterycutout values for each station
SOCAll_df=pd.DataFrame() # to hold all soc for each station
### pv waste ###
pv2wasteAll_df=pd.DataFrame() # to hold all pv waste values for each station
### pv2inv###
pv2invAll_df=pd.DataFrame() # to hold all pv2inv values for each station
### pv2store###
pv2storeAll_df=pd.DataFrame() # to hold all pv2store values for each station

# get the common columns from pv and demand
for stationid in pv_chicago.columns:
    if stationid in demand_all.columns:
        if stationid in rebalance_all.columns:
            pv=pv_chicago[stationid]
            demand=demand_all[stationid]
            rebalance=rebalance_all[stationid]
            E1 = dispatch_max_sc(pv, demand, rebalance,param_tech, return_series=False)
            #results = print_analysis(pv, demand, param_tech, E1)
            BatteryReplaceAll_df[stationid]=E1["BatteryReplace"]
            #ConsumptionAll_df[stationid]=E1['inv2load']
            #BatteryLeftoverAll_df[stationid]=E1['BatteryLeftover']
            BatteryCutoutAll_df[stationid]=E1['BatteryCutout']
            SOCAll_df[stationid]=E1['LevelOfCharge']
            pv2wasteAll_df[stationid] = E1['pv2waste']
            pv2invAll_df[stationid] = E1['pv2inv']
            pv2storeAll_df[stationid] = E1['pv2store']

### save downtime results
BatteryCutoutAll_df.to_csv("results/sys_down_hour.csv")
### pv waste
pv2waste=pv2wasteAll_df.sum(axis=0)
### pv2inv
pv2inv=pv2invAll_df.sum(axis=0)
### pv2store
pv2store=pv2storeAll_df.sum(axis=0)
### pvconsumed0dt
pvconsumedreb=pv2inv*param_tech['InverterEfficiency']+pv2store*param_tech['InverterEfficiency']*param_tech['BatteryEfficiency']
resultsreb=pd.DataFrame()
resultsreb['pv2invreb']=pv2inv
resultsreb['pv2storereb']=pv2store
resultsreb['pv2wastereb']=pv2waste
resultsreb['pvconsumedreb']=pvconsumedreb

resultsreb['pv']=pv_chicago.sum(axis=0)
resultsreb['demand']=demand_all.sum(axis=0)

##SSL
resultsreb['SSLreb']=resultsreb['pvconsumedreb']/resultsreb['demand']*100
## SCL
resultsreb['SCLreb']=resultsreb['pvconsumedreb']/resultsreb['pv']*100
##EPBT CED=3258*n_pv*0.84+1360*n_bat=4096
resultsreb['EPBTreb']=4096/(resultsreb['pvconsumedreb']*3.6/0.3)

###downtime
resultsreb['downtimereb']=BatteryCutoutAll_df.sum(axis=0)/8760*100
##batrep
resultsreb['batrepreb']=BatteryReplaceAll_df.sum(axis=0)

resultsreb.to_csv('results/results_reb.csv')

# ########put results together#########
# df=pd.DataFrame()
# df['batrep']=BatteryReplaceStationTotal
# df['sysavail']=SysAvail
# df['SSL']=ConsumptionEnergyMix['SSL']
# df['eroi']=eroi
# #load inCBD info
# inCBD=pd.read_csv('results/divvystation.csv',index_col=0)
# inCBD.index=inCBD.index.astype(str)
# df_all=pd.merge(df,inCBD,how='inner',left_index=True,right_index=True)
# df_all.to_csv('results/results0220.csv')

############## SOC full ###############################
# #SOCAll_df.to_csv('results/SOCAll_reb.csv')
# fullperc=pd.DataFrame()
# fullperc['full_times']=SOCAll_df[SOCAll_df==1.2].sum(axis=0)


# fullperc['batrep_times']=BatteryReplaceAll_total
# fullperc['socfullonly_times']=fullperc['full_times']-fullperc['batrep_times']
# fullperc['socfullonly_perc']=fullperc['socfullonly_times']/8760*100
#fullperc.to_csv('results/SOC_reb.csv')


############## battery replacement ###############################
#BatteryReplaceAll_df.to_csv('results/BatteryReplaceHour_rebalancing.csv')
# BatteryReplaceMonth=BatteryReplaceAll_df.groupby(BatteryReplaceAll_df.index.month).sum()
# BatteryReplaceMonthMean=BatteryReplaceMonth.mean(axis=1)
#BatteryReplaceMonthMean.to_csv('results/BatteryReplaceMonth_rebalancing.csv')



#BatteryCutoutAll_df.to_csv('results/sys_down_hour.csv')

############## System availabilty ###############################
#cutout percentage
# BatteryCutoutPerc=BatteryCutoutAll_df.sum(axis=0)/BatteryCutoutAll_df.shape[0]*100
# SysAvail=100-BatteryCutoutPerc
# SysAvail=SysAvail.rename('SystemAvailability')
#SysAvail.to_csv('results/SystemAvailability.csv')
#SysAvail.describe()
#boxplot
# plt.boxplot(BatteryCutoutPerc)
# plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom=False,      # ticks along the bottom edge are off
#     top=False,         # ticks along the top edge are off
#     labelbottom=False)
# plt.ylabel('Percentage of Battery Downtime (%)')
# plt.savefig('results/graph/battery_cutout.jpg')

# percentage of cutout during daytime between 6am and 12am
## no big change from results above
# BatteryCutoutDay_df=BatteryCutoutAll_df[BatteryCutoutAll_df.index.hour>=6]
# #cutout percentage
# BatteryCutoutDayPerc=BatteryCutoutDay_df.sum(axis=0)/BatteryCutoutDay_df.shape[0]*100
# BatteryCutoutDayPerc.describe()
# #boxplot
# plt.boxplot(BatteryCutoutDayPerc)
# plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom=False,      # ticks along the bottom edge are off
#     top=False,         # ticks along the top edge are off
#     labelbottom=False)
# plt.ylabel('Percentage of Battery Cutout during Daytime (%)')
# plt.savefig('results/graph/battery_cutout_daytime.jpg')



# ############## SSL/Consumption ###############################
# #ConsumptionAll_df.to_csv('../results/Consumption_hour_1W.csv')
# # annual total by station
# ConsumptionStationTotal=ConsumptionAll_df.sum(axis=0)
# ConsumptionStationTotal=ConsumptionStationTotal.rename("Total")

# ############## BatteryLeftover ###############################
# # annual total by station
# BatteryLeftoverStationTotal=BatteryLeftoverAll_df.sum(axis=0)
# BatteryLeftoverStationTotal=BatteryLeftoverStationTotal.rename("Leftover")
#
# BatteryReplaceStationTotal=BatteryReplaceAll_df.sum(axis=0)
# BatteryReplaceStationTotal=BatteryReplaceStationTotal.rename("ReplacementTimes")

# ############### annual energy mix ###########################
# ## Energy from battery replacement =  batterreplacetimes * BatteryCapacity - batteryleftover
# ## Energy from PV = consumption - energy from battery replacement
# ConsumptionEnergyMix=pd.DataFrame()
# ConsumptionEnergyMix['ConsumptionTotal']=ConsumptionStationTotal
# ConsumptionEnergyMix['BatteryReplacementEnergy']=BatteryReplaceStationTotal*param_tech['BatteryCapacity']*param_tech['InverterEfficiency']*param_tech['BatteryEfficiency']-BatteryLeftoverStationTotal
# ConsumptionEnergyMix['PVEnergy']=ConsumptionEnergyMix['ConsumptionTotal']-ConsumptionEnergyMix['BatteryReplacementEnergy']
# ConsumptionEnergyMix['SSL']=ConsumptionEnergyMix['PVEnergy']/ConsumptionEnergyMix['ConsumptionTotal']*100
# #ConsumptionEnergyMix.to_csv('results/SSL.csv')
# #ConsumptionEnergyMix['SSL'].describe()
#boxplot
#pvpercent_boxplot=ConsumptionEnergyMix['PVperc'].plot.box()
#pvpercent_boxplot.set_ylabel('Percentage of Energy from PV Panel (100%)')
#histogram
#plt.hist(ConsumptionEnergyMix['PVperc'].array, bins='auto')
#ConsumptionEnergyMix[ConsumptionEnergyMix['PVperc']>0.8].count()['PVperc']/ConsumptionEnergyMix.shape[0]

# ###########EROI##################
# CED=4618
# eroi=30*(ConsumptionEnergyMix['PVEnergy']*3.6/0.3)/(CED*0.84)
# eroi_df=pd.DataFrame()
# eroi_df['eroi']=eroi
# #eroi_df.to_csv('results/eroi.csv')


#########################################################
#         monthly result
#########################################################
# battery_replace_hour=pd.read_csv('results/battery_replace_hour_1W.csv',index_col=0)
# times = pd.date_range(start='2019-01-01', end='2019-12-31 23:00:00', freq='H')
# battery_replace_hour=battery_replace_hour.set_index(times)
# battery_replace_month=battery_replace_hour.groupby(battery_replace_hour.index.month).sum()
# battery_replace_month.T.to_csv('results/battery_replace_month_1W.csv')
#BatteryReplaceMonth=pd.read_csv('results/battery_replace_month_1W.csv',index_col=0)
#BatteryReplaceMonthMean=BatteryReplaceMonth.mean(axis=0)
#BatteryReplaceMonthMean.to_csv('results/battery_replace_month_1w_mean.csv')

#########################################################
#         linear relationship
#########################################################
#annual total number of battery replacement by station
# BatteryReplaceAll_df=pd.read_csv('results/battery_replace_hour_1W.csv',index_col=0)
# BatteryReplaceYear=BatteryReplaceAll_df.sum(axis=0)
# BatteryReplaceYear=BatteryReplaceYear.rename('total')
# # join with station capacity
# divvystation=pd.read_csv('data/US_IL_Chicago_Divvy_2021-03-15-11_30_40.csv')
# divvystation=divvystation.set_index('station_id')
# divvystation.index=divvystation.index.map(str) #conver to str to match BatteryReplaceYear
# divvystation_batteryreplaceyear=divvystation.merge(BatteryReplaceYear,how='inner',left_index=True,right_index=True)
# divvystation_batteryreplaceyear.to_csv('results/battery_replace_hour_1W_capacity.csv')
# ConsumptionEnergyMix=pd.read_csv('results/ConsumptionEnergyMix_year_1W.csv',index_col=0)
# ConsumptionEnergyMix.index=ConsumptionEnergyMix.index.map(str)
# divvystation_batteryreplaceyear_energy=divvystation_batteryreplaceyear.merge(ConsumptionEnergyMix,how='inner',left_index=True,right_index=True)
# divvystation_batteryreplaceyear_energy.to_csv('results/divvystation_batteryreplaceyear_energy.csv')
#########################################################
#         PV, Demand and LCA
#########################################################
# pv_chicago_LCA=pv_chicago[BatteryReplaceAll_df.columns]
# #total PV production
# PVStationTotal=pv_chicago_LCA.sum(axis=0)
# #boxplot
# #PVStationTotal.plot.box()
# #histogram
# #plt.hist(PVStationTotal.array, bins='auto')
#
# # how many years to pay back
# #PVStationTotal[PVStationTotal>153.24].count()
# #PVStationTotal[PVStationTotal>76.62].count()
#
# # save years of pay back to PVStationPayback
# # pv values over 152.24, payback in 1 year
# # pv values between 76.62 and 153.24, payback in 2 years
# # otherwise, payback in 3 years or more (year value is 3)
# PVStationPayback = PVStationTotal.to_frame()
# PVStationPayback =PVStationPayback.rename(columns={0:"value"})
# PVStationPayback['year']=0
# PVStationPayback.loc[PVStationPayback['value']>153.24,'year']=1
# PVStationPayback.loc[(PVStationPayback['value']>76.62) & (PVStationPayback['value']<=153.24),'year']=2
# PVStationPayback.loc[PVStationPayback['value']<76.62,'year']=3
# #PVStationPayback.to_csv('../results/pv_station_payback.csv')
#
# #demand
# DemandStationTotal=demand_all.sum(axis=0)

#### SOC #####
# SOC_df=pd.DataFrame() # to hold all SOC values for each station, %
#
# # get the common columns from pv and demand
# for stationid in pv_chicago.columns:
#     if stationid in demand_all.columns:
#         pv=pv_chicago[stationid]
#         demand=demand_all[stationid]
#         E1 = dispatch_max_sc(pv, demand, param_tech, return_series=False)
#         #results = print_analysis(pv, demand, param_tech, E1)
#         SOC_df[stationid]=E1["LevelOfCharge"]/param_tech['BatteryCapacity']
#
# ##count # of 1 in each column/station
# fullperc=pd.DataFrame()
# fullperc['fullperc']=SOC_df[SOC_df==1].sum(axis=0)/8760

#SOC_df.to_csv('results/SOC_all.csv')
#fullperc.to_csv('results/SOC_fullperc.csv')

# plt.boxplot(fullperc['fullperc'])
# plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom=False,      # ticks along the bottom edge are off
#     top=False,         # ticks along the top edge are off
#     labelbottom=False)
# plt.ylabel('percentage with full soc (100%)')