from bright import FRDefaults
fr_params = FRDefaults()

#General Specifications
Quiet = False

#FR Specifications
fr_params.BUt = 27.3969072997					#FR Burnup
fr_params.batches = 1						#Number of FR batches 
fr_params.pnl =  0.60316668602479 				#FR Non-Leakage Probability

#FR Storage
FR_SNF_Storage_Time = "3 y"

#FR Reprocessing
FR_SE_U  = 0.999 
FR_SE_NP = 0.999 
FR_SE_PU = 0.999 
FR_SE_AM = 0.9999 
FR_SE_CM = 0.9999 
FR_SE_CS = 0.999 
FR_SE_SR = 0.999 
