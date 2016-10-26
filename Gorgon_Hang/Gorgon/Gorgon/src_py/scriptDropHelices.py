'''
Script to run drop helices test.
'''
import os

if __name__ == '__main__':

    print '(+++) Drop matches tests start....'

    #
    dataDir = 'H:/Hang/Gorgon/Data/dropHeliesTestDataSet/'
    # ('1ss8_A_1sx4_A','A',19,'2-252','2-252'), ('3145_J_3146_J','J',9,'1-250','1-250'), ('3198_A_3201_A','A',40,'1-926,943-1160','1-926,943-1160'), ('6551_E_6552_A','E',9,'19-349','19-349'), ('6574_K_6575_K','K',12,'1-418','1-418'), ('6574_Q_6575_Q','Q',16,'1-388,390-431','1-388,390-431'), ('6574_Z_6575_Z','Z',31,'99-177,195-993','99-177,195-993'), ('1aon_A_1aon_H','A',17,'2-525','2-525'), ('1oel_A_2c7c_A','A',19,'2-525','2-525'), ('1omp_A_1anf_A','A',16,'1-369','1-369'), ('3tgl_A_4tgl_A','A',8,'5-269','5-269'), ('4ake_A_1ake_A','A',7,'1-214','1-214'), ('9aat_A_1ama_A','A',14,'3-410','3-410'),('3izh_C_3j3x_I','C',18,'1031-1538','11-518'), ('2c7c_M_2c7c_A','M',20,'3-524','3-524')
    sampleList = [('6551_E_6552_A','E',9,'19-349','19-349')];

    for sample in sampleList:
        calphaAtomFile = dataDir+sample[0]+'/in.pdb'
        # Well, we have two possible file format
        helixFile = dataDir+sample[0]+'/in.wrl'
        if not os.path.isfile(helixFile):
            helixFile = dataDir+sample[0]+'/in.vrml'
            
        skeletonFile = dataDir+sample[0]+'/in.off'
        correspondenceFile = dataDir+sample[0]+'/in.cor'
        chainID  = sample[1]
        matchNum = sample[2]
        # We keep dropping until only two helices left
        for i in range(1, matchNum-1):
            dropHelicesNum = str(i)
            outputFile = dataDir + sample[0] + '/output/' + str(i) + '.pdb'
            os.system('python gorgon.pyw '+calphaAtomFile+' '+chainID+' '+helixFile+' '+skeletonFile+' '+correspondenceFile+' '+dropHelicesNum+' '+outputFile)   
            
    '''
        calphaAtomFile = 'H:/Hang/Gorgon/Data/Derek_Real_Data_copy_2/2c7c_M_To_2c7c_A/2c7c_M.pdb'
        chainID = 'M'
        helixFile = 'H:\Hang\Gorgon\Data\Derek_Real_Data_copy_2\\2c7c_M_To_2c7c_A\\2c7c_A_hunter.wrl'
        skeletonFile = 'H:\Hang\Gorgon\Data\Derek_Real_Data_copy_2\\2c7c_M_To_2c7c_A\\2c7c_A.off'
        correspondenceFile = 'H:\Hang\Gorgon\Data\Derek_Real_Data_copy_2\\2c7c_M_To_2c7c_A\\2c7c_M_2c7c_A.cor'
        outputName = 'H:\Hang\Gorgon\Data\Derek_Real_Data_copy_2\\2c7c_M_To_2c7c_A\\testOutput.pdb'
        matchNum = 2
        for i in range(1, matchNum+1):
            dropHelicesNum = str(i)
            #outputName = outputDir + 
            os.system('python gorgon.pyw '+calphaAtomFile+' '+chainID+' '+helixFile+' '+skeletonFile+' '+correspondenceFile+' '+dropHelicesNum+' '+outputName)
    '''

    print '(+++) Done dropping matches tests!'
        