import pandas as pd
import numpy as np
import os

##### final/intermediate results saved ###
#customer_count.to_csv('../results/customer_count2019.csv')
#bike_count.to_csv('../results/bike_count2019.csv')

#2019 missing station filled with GBFS data first and then 2021 data

######### calculate demand ###########
# demand = # bikes docked * 1W + # trips with casual users * 30W * 4min/60 (or max 30W)


######################################################
################### divvy trip data #################
######################################################
# customer trip start
############ process Q1 ###############
#load trip data Q1
divvyQ1 = pd.read_csv('../data/divvytrip/Divvy_Trips_2019_Q1.csv')
#convert start_time amd end_time to datetime
divvyQ1['start_time_dt']=pd.to_datetime(divvyQ1['start_time']).dt.floor('H')

############ process Q2 ###############
#load trip data Q2
divvyQ2 = pd.read_csv('../data/divvytrip/Divvy_Trips_2019_Q2.csv')
#convert start_time amd end_time to datetime
divvyQ2['start_time_dt']=pd.to_datetime(divvyQ2['01 - Rental Details Local Start Time']).dt.floor('H')

############ process Q3 ###############
#load trip data Q3
divvyQ3 = pd.read_csv('../data/divvytrip/Divvy_Trips_2019_Q3.csv')
#convert start_time amd end_time to datetime
divvyQ3['start_time_dt']=pd.to_datetime(divvyQ3['start_time']).dt.floor('H')

############ process Q4 ###############
#load trip data Q4
divvyQ4 = pd.read_csv('../data/divvytrip/Divvy_Trips_2019_Q4.csv')
#convert start_time amd end_time to datetime
divvyQ4['start_time_dt']=pd.to_datetime(divvyQ4['start_time']).dt.floor('H')

#concat 4 df
divvy_df=pd.concat([divvyQ1[['start_time_dt','from_station_id','trip_id','usertype']],
                    divvyQ2[['start_time_dt','03 - Rental Start Station ID','01 - Rental Details Rental ID','User Type']].rename(columns={'03 - Rental Start Station ID':'from_station_id','01 - Rental Details Rental ID':'trip_id','User Type':'usertype'}),
                    divvyQ3[['start_time_dt','from_station_id','trip_id','usertype']],
                    divvyQ4[['start_time_dt','from_station_id','trip_id','usertype']]],ignore_index=True)

#time index
#keep only 2019 dates
times_df = pd.DataFrame(pd.date_range(start='2019-01-01', end='2019-12-31 23:00:00', freq='H'),columns=['date'])

########## user type ###########
## keep only customer ## not subscriber
customer_df=divvy_df.loc[divvy_df['usertype']=='Customer']
customer_trip=customer_df.groupby(['start_time_dt','from_station_id']).count()
customer_trip=customer_trip.reset_index()
customer_pivot=customer_trip.pivot(index="start_time_dt",columns="from_station_id",values="trip_id")
customer_count=times_df.merge(customer_pivot,left_on='date',right_index=True,how='left')
customer_count=customer_count.set_index('date')
customer_count=customer_count.fillna(0) #fill nan as 0
#customer_trip.iloc[customer_trip['trip_id'].idxmax()] # max 121, station 35, 2019-8-3 16:00:00
#customer_trip['trip_id'].describe() #mean 2.43 75% 2 max 121

# customer_count is the hourly count of # customers by station
#customer_count.to_csv('../results/customer_count2019.csv')

######################################################
##################divvy station to calculate # bikes docked by hour#####################
######################################################
divvy_station_2019_raw = pd.read_csv('../data/divvystation/Divvy_Bicycle_Stations2019.csv')
divvy_station_2019_raw['hour'] = pd.to_datetime(divvy_station_2019_raw['Timestamp']).dt.floor('H')
divvy_station_2019 = divvy_station_2019_raw.groupby(['ID', 'hour']).mean().reset_index()
divvy_Station_2019_pivot = divvy_station_2019.pivot(index="hour", columns="ID", values="Available Bikes")
divvy_Station_2019_pivot=divvy_Station_2019_pivot.fillna(0)  #fill nan as 0

# find missing timestamp
divvy_station_2019_na=times_df.merge(divvy_Station_2019_pivot,left_on='date',right_index=True,how='left')[times_df.merge(divvy_Station_2019_pivot,left_on='date',right_index=True,how='left').isnull().any(axis=1)]
divvy_station_2019_na=divvy_station_2019_na.set_index('date')  #datetime as index
divvy_station_2019_na=divvy_station_2019_na.drop(columns=divvy_station_2019_na.columns)  #keep only the datetime
#missing dates list (as string)
divvy_station_2019_na_datestr=pd.to_datetime(divvy_station_2019_na.index.date,format='%Y-%m-%d').astype(str).unique()
divvy_station_2019_na_datestr_gbfsfolder='Chicago_'+divvy_station_2019_na_datestr
divvy_station_2019_na_datestr_gbfsfolder=divvy_station_2019_na_datestr_gbfsfolder.to_series().reset_index()
#divvy_station_2019_na_datestr_gbfsfolder.to_csv('data/divvystation/divvy2019_na_checkGBFS.csv')
#data found in gbfs file
divvy_station_2019_na_gbfs_raw=pd.read_csv('../data/divvystation/divvy2019_na_gbfs.csv')
#field Download Time is the time
#field availableBikes and num_bikes_available same, copy availableBikes to num_bikes_available
#station_id and id are the same, copy id to station_id
# so the fields used are Download Time, num_bikes_available and station_id
divvy_station_2019_na_gbfs_raw.loc[divvy_station_2019_na_gbfs_raw['station_id'].isnull(),'station_id']=divvy_station_2019_na_gbfs_raw['id']
divvy_station_2019_na_gbfs_raw.loc[divvy_station_2019_na_gbfs_raw['num_bikes_available'].isnull(),'num_bikes_available']=divvy_station_2019_na_gbfs_raw['availableBikes']
divvy_station_2019_na_gbfs_raw['station_id']=divvy_station_2019_na_gbfs_raw['station_id'].astype(int)
divvy_station_2019_na_gbfs_raw['Hour']=pd.to_datetime(divvy_station_2019_na_gbfs_raw['Download Time']).dt.floor('H')
divvy_station_2019_na_gbfs = divvy_station_2019_na_gbfs_raw.groupby(['station_id', 'Hour']).mean().reset_index()
divvy_station_2019_na_gbfs_pivot=divvy_station_2019_na_gbfs.pivot(index="Hour", columns="station_id", values="num_bikes_available")
#drop columns in gbfs but not in 2019
common_cols = list(set(divvy_Station_2019_pivot.columns).intersection(divvy_station_2019_na_gbfs_pivot.columns)) #360, 363, 397 not in 2020 but in 2019
divvy_station_2019_na_gbfs_pivot = divvy_station_2019_na_gbfs_pivot[common_cols]
divvy_station_2019_na_gbfs_pivot=divvy_station_2019_na_gbfs_pivot.fillna(0)
#join divvy_station_2019_na_gbfs_pivot to na
divvy_station_2019_na_gbfs_filled=divvy_station_2019_na.merge(divvy_station_2019_na_gbfs_pivot,left_index=True,right_index=True,how="left")

#fill the remaining na with 2021 data
divvy_station_2019_na_gbfs_filled_na=divvy_station_2019_na_gbfs_filled.loc[divvy_station_2019_na_gbfs_filled[2].isnull()].index
divvy_station_2021_raw = pd.read_csv('../data/divvystation/Divvy_Bicycle_Stations20210101_1213.csv')
divvy_station_2021_raw['hour'] = pd.to_datetime(divvy_station_2021_raw['Timestamp']).dt.floor('H')
divvy_station_2021 = divvy_station_2021_raw.groupby(['ID', 'hour']).mean().reset_index()
divvy_Station_2021_pivot = divvy_station_2021.pivot(index="hour", columns="ID", values="Available Bikes")
divvy_Station_2021_pivot=divvy_Station_2021_pivot.fillna(0)  #fill nan as 0

#2021 times
times_df_2021 = pd.DataFrame(pd.date_range(start='2021-01-01', end='2021-12-31 23:00:00', freq='H'),columns=['date'])
#keep only 2021 dates
divvy_station_2021_pivot=times_df_2021.merge(divvy_Station_2021_pivot,left_on='date',right_index=True,how='left')
#deduct a year for 2020 data to join with 2019
divvy_station_2021_pivot['date_offset']=divvy_station_2021_pivot['date']+pd.offsets.DateOffset(years=-2)
divvy_station_2021_pivot = divvy_station_2021_pivot.set_index('date_offset')
divvy_station_2021_pivot = divvy_station_2021_pivot.drop(columns='date')
#drop columns in 2021 but not in 2019
common_cols2021=list(set(divvy_station_2021_pivot.columns).intersection(common_cols))
divvy_station_2021_pivot_offset = divvy_station_2021_pivot[common_cols2021]
#join 2021 offset to 2019_na_gbfs_na
divvy_station_2019_na_gbfs_filled_na_df=pd.DataFrame(index=divvy_station_2019_na_gbfs_filled_na)
divvy_station_2019_na_gbfs_filled_na_filled=divvy_station_2019_na_gbfs_filled_na_df.merge(divvy_station_2021_pivot_offset,left_index=True,right_index=True,how="left")

#drop columns
divvy_Station_2019_pivot=divvy_Station_2019_pivot[common_cols2021]
divvy_station_2019_na_gbfs_filled_notna=divvy_station_2019_na_gbfs_filled[divvy_station_2019_na_gbfs_filled[2].notnull()]
divvy_station_2019_na_gbfs_filled_notna=divvy_station_2019_na_gbfs_filled_notna[common_cols2021]
#append to divvy 2019 data
bike_count = divvy_Station_2019_pivot.append(divvy_station_2019_na_gbfs_filled_notna).append(divvy_station_2019_na_gbfs_filled_na_filled)
bike_count=bike_count.fillna(0)
#bike_count.to_csv('../results/bike_count2019.csv')

######################################################
################## calculate demand #####################
######################################################
# demand = # bikes docked * 1W + # trips with casual users * 30W * 4min/60 (or max 30W) that is *2
#unit kw
demand=pd.DataFrame()
for stationid in customer_count.columns:
    if stationid in bike_count.columns:
        demand[stationid] = (10+bike_count[stationid] + np.minimum(customer_count[stationid]*2,30))*0.001

#demand.to_csv('../results/demand2019.csv')
##annual demand
demand = pd.read_csv('results/demand2019_1W.csv',index_col=0)
demand_year=demand.sum(axis=0)
demand_year.hist()
plt.xlabel('Annual demand (kWh)')
plt.ylabel('Number of stations')
plt.show()

#demand.mean().mean() # 0.0072kW on average per hour

#daily demand
times = pd.date_range(start='2019-01-01', end='2019-12-31 23:00:00', freq='H')
demand=demand.set_axis(times)
demand_day = demand.resample('D').sum()
#demand_day.mean().mean() # 0.414kWH
#demand.resample('D').sum().mean(axis=1).hist()

#monthly demand
demand_month = demand.resample('M').sum()
#demand_month.mean(axis=1)

####### demand breakdown ############
times = pd.date_range(start='2019-01-01', end='2019-12-31 23:00:00', freq='H')
customer_count=pd.read_csv('results/customer_count2019.csv', index_col=0)
bike_count=pd.read_csv('results/bike_count2019.csv', index_col=0)
demand_customer=pd.DataFrame()
for stationid in customer_count.columns:
    if stationid in bike_count.columns:
        demand_customer[stationid] = np.minimum(customer_count[stationid]*2,30)*0.001
#daily demand from customers/kiosk use
demand_customer=demand_customer.set_axis(times)
demand_day_customer = demand_customer.resample('D').sum()
demand_day_customer_mean=demand_day_customer.mean(axis=1)
demand_customer_year=demand_customer.sum(axis=0)

demand_bike=pd.DataFrame()
for stationid in customer_count.columns:
    if stationid in bike_count.columns:
        demand_bike[stationid] = (10+bike_count[stationid])*0.001
#daily demand from bike
demand_bike=demand_bike.set_axis(times)
demand_day_bike = demand_bike.resample('D').sum()
demand_day_bike_mean=demand_day_bike.mean(axis=1)
demand_bike_year=demand_bike.sum(axis=0)

demand_year_breakdown=pd.DataFrame()
demand_year_breakdown['Kiosk']=demand_customer_year
demand_year_breakdown['Bike Dock']=demand_bike_year

plt.hist(demand_year_breakdown, 20, density=False, histtype='bar', stacked=True)
plt.xlabel('Annual demand (kWh)')
plt.ylabel('Number of stations')
plt.show()