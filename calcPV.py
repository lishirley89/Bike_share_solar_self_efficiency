# calculate pv profile for chicago all stations
import pandas as pd
import numpy as np

#load radiation for chicago.
# dataframe row 0-9, 1/15; 10-20, 2/15 ; 21-32, 3/15; 33-46, 4/15
# 47-61, 5/15; 62-76, 6/15; 77-91, 7/15; 92-105, 8/15;
# 106-118, 9/13; 119-129, 10/15; 130-139, 11/15; 140-148, 12/15
rad_chicago = pd.read_csv('../data/rad_chicago_year.csv')
#calcualte pv from radiation, assume 15% efficiency and 86% performance ration, kWh
# total area of pv panel, 0.84m2
pv_chicago_raw = rad_chicago.multiply(0.001*0.15*0.86*0.84)
# copy pv_chicago_raw and insert 0 or the pv values
pv_chicago = pd.DataFrame(columns=pv_chicago_raw.columns)

#jan
pv01=pd.DataFrame(0, index=np.arange(7), columns=pv_chicago.columns)
pv01=pv01.append(pv_chicago_raw.loc[0:9])
pv01=pv01.append(pd.DataFrame(0, index=np.arange(7), columns=pv_chicago.columns))
for i in range(31):
    pv_chicago=pv_chicago.append(pv01)

#feb
pv02=pd.DataFrame(0, index=np.arange(6), columns=pv_chicago.columns)
pv02=pv02.append(pv_chicago_raw.loc[10:20])
pv02=pv02.append(pd.DataFrame(0, index=np.arange(7), columns=pv_chicago.columns))
for i in range(28):
    pv_chicago=pv_chicago.append(pv02)

#mar
pv03=pd.DataFrame(0, index=np.arange(6), columns=pv_chicago.columns)
pv03=pv03.append(pv_chicago_raw.loc[21:32])
pv03=pv03.append(pd.DataFrame(0, index=np.arange(6), columns=pv_chicago.columns))
for i in range(31):
    pv_chicago=pv_chicago.append(pv03)

#apr
pv04=pd.DataFrame(0, index=np.arange(6), columns=pv_chicago.columns)
pv04=pv04.append(pv_chicago_raw.loc[33:46])
pv04=pv04.append(pd.DataFrame(0, index=np.arange(4), columns=pv_chicago.columns))
for i in range(30):
    pv_chicago=pv_chicago.append(pv04)

#may
pv05=pd.DataFrame(0, index=np.arange(5), columns=pv_chicago.columns)
pv05=pv05.append(pv_chicago_raw.loc[47:61])
pv05=pv05.append(pd.DataFrame(0, index=np.arange(4), columns=pv_chicago.columns))
for i in range(31):
    pv_chicago=pv_chicago.append(pv05)

#jun
pv06=pd.DataFrame(0, index=np.arange(5), columns=pv_chicago.columns)
pv06=pv06.append(pv_chicago_raw.loc[62:76])
pv06=pv06.append(pd.DataFrame(0, index=np.arange(4), columns=pv_chicago.columns))
for i in range(30):
    pv_chicago=pv_chicago.append(pv06)

#jul
pv07=pd.DataFrame(0, index=np.arange(5), columns=pv_chicago.columns)
pv07=pv07.append(pv_chicago_raw.loc[77:91])
pv07=pv07.append(pd.DataFrame(0, index=np.arange(4), columns=pv_chicago.columns))
for i in range(31):
    pv_chicago=pv_chicago.append(pv07)

#aug
pv08=pd.DataFrame(0, index=np.arange(5), columns=pv_chicago.columns)
pv08=pv08.append(pv_chicago_raw.loc[92:105])
pv08=pv08.append(pd.DataFrame(0, index=np.arange(5), columns=pv_chicago.columns))
for i in range(31):
    pv_chicago=pv_chicago.append(pv08)

#sep
pv09=pd.DataFrame(0, index=np.arange(7), columns=pv_chicago.columns)
pv09=pv09.append(pv_chicago_raw.loc[106:118])
pv09=pv09.append(pd.DataFrame(0, index=np.arange(4), columns=pv_chicago.columns))
for i in range(30):
    pv_chicago=pv_chicago.append(pv09)

#oct
pv10=pd.DataFrame(0, index=np.arange(7), columns=pv_chicago.columns)
pv10=pv10.append(pv_chicago_raw.loc[119:129])
pv10=pv10.append(pd.DataFrame(0, index=np.arange(6), columns=pv_chicago.columns))
for i in range(31):
    pv_chicago=pv_chicago.append(pv10)

#nov
pv11=pd.DataFrame(0, index=np.arange(7), columns=pv_chicago.columns)
pv11=pv11.append(pv_chicago_raw.loc[130:139])
pv11=pv11.append(pd.DataFrame(0, index=np.arange(7), columns=pv_chicago.columns))
for i in range(30):
    pv_chicago=pv_chicago.append(pv11)

#dec
pv12=pd.DataFrame(0, index=np.arange(8), columns=pv_chicago.columns)
pv12=pv12.append(pv_chicago_raw.loc[140:148])
pv12=pv12.append(pd.DataFrame(0, index=np.arange(7), columns=pv_chicago.columns))
for i in range(31):
    pv_chicago=pv_chicago.append(pv12)

#set index
times = pd.date_range(start='2019-01-01', end='2019-12-31 23:00:00', freq='H')
pv_chicago.set_index(times,inplace=True)
pv_chicago.to_csv('../data/pv_chicago_year.csv')

# get annual solar radiation distribution
#pv_chicago_hour=pd.read_csv('data/pv_chicago_year.csv', index_col=0)
pv_chicago_year=pv_chicago.sum(axis=0)
pv_chicago_year_boxplot=pv_chicago_year.plot.box()
pv_chicago_year_boxplot.set_ylabel('Annual PV Production (kWh/m2)')
#pv_chicago_year.to_csv('../results/pv_chicago_yeartotal.csv')
rad_chicago_year=pv_chicago_year/0.108  # 0.15*0.86*0.84=0.108, unit kwh/m2
#rad_chicago_year_boxplot=rad_chicago_year.plot.box()
#rad_chicago_year_boxplot.set_ylabel('Annual Solar Radiation (kWh/m2)')
