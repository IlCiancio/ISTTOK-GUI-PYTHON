"""
Created on Mon Feb  4 17:58:27 2019

@author: Il_Ciancio
"""
#SOME FUNCTION ON THE DATABASE

def getValueFromChannel(Client,Channel,ShotNumber):
	print("get Value from: " + Channel)
	import matplotlib.pyplot as plt
	from sdas.core.SDAStime import Date, Time, TimeStamp
	import numpy as np
	structArray = Client.getData(Channel, '0x0000', ShotNumber)
	struct = structArray[0]
	Value = struct.getData()
	#TIME
	tstart = struct.getTStart()
	tend = struct.getTEnd()
    #Calculate the time between samples
	tbs = ( tend.getTimeInMicros() - tstart.getTimeInMicros() )/(len(Value)*1.0)
	#Get the events  associated with this data
	events = struct.get('events')
	tevent = TimeStamp(tstamp=events[0].get('tstamp')) #The delay of the start time relative to the event time
	delay = tstart.getTimeInMicros() - tevent.getTimeInMicros()
	#Finally create the time array
	times = np.linspace(delay,delay+tbs*(len(Value)-1),len(Value))
	plt.plot(times, Value)
	plt.xlabel('time Micros')
	plt.ylabel(Channel)
	plt.grid(True)
	plt.show()
	return Value, times
	
def plotValueFromChannel(Channel,ShotNumber):
    print(Channel)
    from sdas.core.client.SDASClient import SDASClient
    import sys

    def StartSdas():	
        host = 'baco.ipfn.ist.utl.pt'
        port = 8888
        client = SDASClient(host,port)	
        return client
    print("Import data relative to the ShotNumber:" + str(ShotNumber))	
    client = StartSdas()
    data = getValueFromChannel(client,Channel,ShotNumber)
    return data
    


def Bmagnmirnv( Z_filament,R_filament,I_filament,r_mirnv,z_mirnv):
#-------------------------------------------------------------------------#

#---Filament is in the R-Y plane and Magnetic Field is Evaluated -------------#

#-------------at the Mirnv coordinate in the R-Z plane------------------------#

#-------------------------------------------------------------------------#
	import numpy as np
	from numpy import pi
	Zc = Z_filament
	I = I_filament
	turns = 1															 # I just have one filament of one turn
	N = 100  								   							 # No of grids in the coil ( X-Y plane)
	u0 = 4*pi*0.001  													 # [microWb/(A cm)]
	phi = np.linspace(0, 2*pi, N) 									 	 # For describing a circle (coil)
	Rc = R_filament * np.cos(phi) 										 #R-coordinates of the filament
	Yc = R_filament * np.sin(phi) 										 #Y-coordinates of the filament
	
	#PYTHON RETURN ERROR IF I DON'T PREALLOCATE A VECTOR
	#Lets obtain the position vectors from dl to the mirnov 
	#mirnov is localized in the plane (y=0)

	RR = np.zeros(N)
	Rz = np.zeros(N)
	Ry = np.zeros(N)
	dlR = np.zeros(N)
	dly = np.zeros(N)
	for i in range(N-1):
		RR[i] = r_mirnv - 0.5* (Rc[i]+Rc[i+1]) 
		Rz[i] = z_mirnv - Zc
		Ry[i] = -0.5 * (Yc[i]+Yc[i+1])
		dlR[i]= Rc[i+1]-Rc[i]
		dly[i]=Yc[i+1]-Yc[i]
		
	RR[-1] = r_mirnv - 0.5*(Rc[-1]+Rc[0])
	Rz[-1] = z_mirnv - Zc
	Ry[-1] = -0.5 * (Rc[-1] + Rc[0])
	dlR[-1] = -Rc[-1] + Rc[0]
	dly[-1] = -Yc[-1]+Yc[0]
	
	#dl x r
	Rcross = np.multiply(-dly,Rz) #or -dly * Rz
	Ycross = np.multiply(dlR,Rz)
	Zcross = np.multiply(dly,RR) - np.multiply(dlR,Ry)
	R = np.sqrt(np.square(RR) + np.square(Rz) + np.square(Ry)) #OR 	#R = sqrt(RR**2 + Rz**2 + Ry**2) #R = np.sqrt(RR**2 + Rz**2 + Ry**2)
	
	#dB=m0/4pi (Idl x r)/r^2 
	BR1 = np.multiply(np.divide(I*u0, 4*pi*R**3), Rcross)
	Bz1 = np.multiply(np.divide(I*u0, 4*pi*R**3), Zcross)
	By1 = np.multiply(np.divide(I*u0, 4*pi*R**3), Ycross)
	
	#Initialize sum magnetic field to be zero first
	BR = 0
	Bz = 0
	By = 0
	BR = BR + np.sum(BR1)
	Bz = Bz + np.sum(Bz1)
	By = By + np.sum(By1)
	
	BR = BR * turns
	By = By * turns
	Bz = Bz * turns #units=[uWb / cm^2]
	
	vector = [Z_filament - z_mirnv, R_filament - r_mirnv]  #Vector from center of chamber to mirnov center
	unit_vec = np.divide(vector, np.linalg.norm(vector) )  # Unit vector
	norm_vec = [unit_vec[1], -unit_vec[0] ] 			   # Normal vector, coil direction
	Bmirn = np.absolute( BR * unit_vec[1] + Bz * unit_vec[0])
	Bmirn = np.dot([Bz,BR], norm_vec)
	Bmirn = 0.01 * Bmirn 									 #fator de 0.01 pra ter [T] 
	return Bmirn


	
def getDataFromDatabase(ShotNumber):
	from sdas.core.client.SDASClient import SDASClient
	from sdas.core.SDAStime import Date, Time, TimeStamp
	import numpy as np
	import matplotlib.pyplot as plt
	import sys
    
	def StartSdas():
		host = 'baco.ipfn.ist.utl.pt'
		port = 8888
		client = SDASClient(host,port)	
		return client

	#CHANGE SHOT NUMBER HERE
	print("Import data relative to the ShotNumber:" + str(ShotNumber))
	shotnr = ShotNumber	
	client = StartSdas()

	#numerical integrated mirnov coils signals
	mirnv =	[	'MARTE_NODE_IVO3.DataCollection.Channel_129',
				'MARTE_NODE_IVO3.DataCollection.Channel_130',
				'MARTE_NODE_IVO3.DataCollection.Channel_131',
				'MARTE_NODE_IVO3.DataCollection.Channel_132',
				'MARTE_NODE_IVO3.DataCollection.Channel_133',
				'MARTE_NODE_IVO3.DataCollection.Channel_134',
				'MARTE_NODE_IVO3.DataCollection.Channel_135',
				'MARTE_NODE_IVO3.DataCollection.Channel_136',
				'MARTE_NODE_IVO3.DataCollection.Channel_137',
				'MARTE_NODE_IVO3.DataCollection.Channel_138',
				'MARTE_NODE_IVO3.DataCollection.Channel_139',
				'MARTE_NODE_IVO3.DataCollection.Channel_140']

	# mirnv_raw=['MARTE_NODE_IVO3.DataCollection.Channel_145',
	# 'MARTE_NODE_IVO3.DataCollection.Channel_146',
	# 'MARTE_NODE_IVO3.DataCollection.Channel_147',
	# 'MARTE_NODE_IVO3.DataCollection.Channel_148',
	# 'MARTE_NODE_IVO3.DataCollection.Channel_149',
	# 'MARTE_NODE_IVO3.DataCollection.Channel_150',
	# 'MARTE_NODE_IVO3.DataCollection.Channel_151',
	# 'MARTE_NODE_IVO3.DataCollection.Channel_152',
	# 'MARTE_NODE_IVO3.DataCollection.Channel_153',
	# 'MARTE_NODE_IVO3.DataCollection.Channel_154',
	# 'MARTE_NODE_IVO3.DataCollection.Channel_155',
	# 'MARTE_NODE_IVO3.DataCollection.Channel_156']

	#numerical integrated mirnov coils signals with offset correction
	mirnv_corr = [	'MARTE_NODE_IVO3.DataCollection.Channel_166',
					'MARTE_NODE_IVO3.DataCollection.Channel_167',
					'MARTE_NODE_IVO3.DataCollection.Channel_168',
					'MARTE_NODE_IVO3.DataCollection.Channel_169',
					'MARTE_NODE_IVO3.DataCollection.Channel_170',
					'MARTE_NODE_IVO3.DataCollection.Channel_171',
					'MARTE_NODE_IVO3.DataCollection.Channel_172',
					'MARTE_NODE_IVO3.DataCollection.Channel_173',
					'MARTE_NODE_IVO3.DataCollection.Channel_174',
					'MARTE_NODE_IVO3.DataCollection.Channel_175',
					'MARTE_NODE_IVO3.DataCollection.Channel_176',
					'MARTE_NODE_IVO3.DataCollection.Channel_177']

	#Calculated external flux on each minov coil through the SS modele
	ext_flux = ['MARTE_NODE_IVO3.DataCollection.Channel_214',
				'MARTE_NODE_IVO3.DataCollection.Channel_215',
				'MARTE_NODE_IVO3.DataCollection.Channel_216',
				'MARTE_NODE_IVO3.DataCollection.Channel_217',
				'MARTE_NODE_IVO3.DataCollection.Channel_218',
				'MARTE_NODE_IVO3.DataCollection.Channel_219',
				'MARTE_NODE_IVO3.DataCollection.Channel_220',
				'MARTE_NODE_IVO3.DataCollection.Channel_221',
				'MARTE_NODE_IVO3.DataCollection.Channel_222',
				'MARTE_NODE_IVO3.DataCollection.Channel_223',
				'MARTE_NODE_IVO3.DataCollection.Channel_224',
				'MARTE_NODE_IVO3.DataCollection.Channel_225']

	#Minov Coils signals with the effect from the external fluxes subtracted		
	mirnv_corr_flux = [	'MARTE_NODE_IVO3.DataCollection.Channel_202',
						'MARTE_NODE_IVO3.DataCollection.Channel_203',
						'MARTE_NODE_IVO3.DataCollection.Channel_204',
						'MARTE_NODE_IVO3.DataCollection.Channel_205',
						'MARTE_NODE_IVO3.DataCollection.Channel_206',
						'MARTE_NODE_IVO3.DataCollection.Channel_207',
						'MARTE_NODE_IVO3.DataCollection.Channel_208',
						'MARTE_NODE_IVO3.DataCollection.Channel_209',
						'MARTE_NODE_IVO3.DataCollection.Channel_210',
						'MARTE_NODE_IVO3.DataCollection.Channel_211',
						'MARTE_NODE_IVO3.DataCollection.Channel_212',
						'MARTE_NODE_IVO3.DataCollection.Channel_213']

	#Measured currents applied by the Primary, Horizontal and Vertical PowerSupply						
	prim = 'MARTE_NODE_IVO3.DataCollection.Channel_093'
	hor = 'MARTE_NODE_IVO3.DataCollection.Channel_091'
	vert = 'MARTE_NODE_IVO3.DataCollection.Channel_092'
    
	#plasma current measured by the Rogowski coil
	Ip_rog = 'MARTE_NODE_IVO3.DataCollection.Channel_088'
	Ip_rog_value = getValueFromChannel(client,Ip_rog,shotnr)
	
    #chopper =' MARTE_NODE_IVO3.DataCollection.Channel_141'

    #Plasma current reconstructed by the mirnov coils without the correction from external fluxes
	Ip_magn = 'MARTE_NODE_IVO3.DataCollection.Channel_085'
	Ip_magn_value = getValueFromChannel(client,Ip_magn,shotnr)

    #Plasma current reconstructed by the mirnov coils without the correction from external fluxes
	Ip_magn_corr = 'MARTE_NODE_IVO3.DataCollection.Channel_228'
	Ip_magn_corr_value = getValueFromChannel(client,Ip_magn_corr,shotnr)


	#SAVES MIRNOV DATA IN A MATRIX
	coilNr=0
	data=[]
	for coil in mirnv_corr:
		coilNr+=1
		structArray=client.getData(coil,'0x0000', shotnr)
		struct=structArray[0]
		data.append(struct.getData())

	#TIME
	tstart = struct.getTStart()
	tend = struct.getTEnd()
    
	#Calculate the time between samples
	tbs = (tend.getTimeInMicros() - tstart.getTimeInMicros())/(len(data[coilNr-1])*1.0)
	
    #Get the events  associated with this data
	events = struct.get('events')
	tevent = TimeStamp(tstamp=events[0].get('tstamp'))
	
    #The delay of the start time relative to the event time
	delay = tstart.getTimeInMicros() - tevent.getTimeInMicros()
	
    #Finally create the time array
	times = np.linspace(delay,delay+tbs*(len(data[coilNr-1])-1),len(data[coilNr-1]))
    #plt.plot(times, data[coilNr-1]);

    #PLOTS ALL DATA FROM MIRNOVS
	coilNr=0
	for coil in data:
		coilNr+=1
		#plt.title('Coil' + str(coilNr))
		ax = plt.subplot(4, 3, coilNr)
		ax.set_title('coil'+str(coilNr))
		plt.plot(times, coil)
		plt.grid(True)
		
	plt.show()
	
	return [Ip_magn_value, Ip_magn_corr_value, Ip_rog_value]
	
	
def getDataForGUI(ShotNumber):
	from sdas.core.client.SDASClient import SDASClient
	from sdas.core.SDAStime import Date, Time, TimeStamp
	import numpy as np
	import matplotlib.pyplot as plt
	import sys
    
	def StartSdas():
		host = 'baco.ipfn.ist.utl.pt'
		port = 8888
		client = SDASClient(host,port)	
		return client

	#CHANGE SHOT NUMBER HERE
	print("Import data relative to the ShotNumber:" + str(ShotNumber))
	shotnr = ShotNumber	
	client = StartSdas()


	#POSITION OF CENTROID					
	R_fromMirnov = 'MARTE_NODE_IVO3.DataCollection.Channel_083'
	R_fromMirnov_value, t1 = getValueFromChannel(client,R_fromMirnov,shotnr)
	
	z_fromMirnov = 'MARTE_NODE_IVO3.DataCollection.Channel_084'
	z_fromMirnov_value, t2 = getValueFromChannel(client,z_fromMirnov,shotnr)
	
	R_fromProbes = 'MARTE_NODE_IVO3.DataCollection.Channel_081'
	R_fromProbes_value, t3 = getValueFromChannel(client,R_fromProbes,shotnr)
	
	z_fromProbes = 'MARTE_NODE_IVO3.DataCollection.Channel_082'
	z_fromProbes_value, t4 = getValueFromChannel(client,z_fromProbes,shotnr)
        
	
	return R_fromMirnov_value, z_fromMirnov_value, R_fromProbes_value, z_fromProbes_value, t1