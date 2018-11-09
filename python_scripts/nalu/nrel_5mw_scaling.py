import os
import numpy as np
import matplotlib.pyplot as plt
try:
    import python_scripts.nalu.io as nalu
except ImportError:
    raise ImportError('Download https://github.com/lawsonro3/python_scripts/blob/master/python_scripts/nalu/nalu_functions.py')

if __name__ == '__main__':

    root_dir = '/Users/mlawson/GoogleDrive/Work/NREL/Projects/HFM-ECP/nrel_5mw/results/cori_data/'

    if os.path.isdir(root_dir) is False:
        raise Exception('root_dir does not exist')

    ####################################
    # Load gC data
    ####################################
    file_gC_13 = root_dir+'gCoarse.13/nrel_5mw_gCoarse.log'
    th_gC_13,t_gC_13 = nalu.read_log(file_gC_13)
    t_gC_13_avg = np.mean(t_gC_13[375:425,:],axis=0)

    file_gC_26 = root_dir+'gCoarse.26/nrel_5mw_gCoarse.log'
    th_gC_26,t_gC_26 = nalu.read_log(file_gC_26)
    t_gC_26_avg = np.mean(t_gC_26[300:350,:],axis=0)

    file_gC_52 = root_dir+'gCoarse.52/nrel_5mw_gCoarse.log'
    th_gC_52,t_gC_52 = nalu.read_log(file_gC_52)
    t_gC_52_avg = np.mean(t_gC_52[500:550,:],axis=0)

    file_gC_104 = root_dir+'gCoarse.104/nrel_5mw_gCoarse.log'
    th_gC_104,t_gC_104 = nalu.read_log(file_gC_104)
    t_gC_104_avg = np.mean(t_gC_104[200:250,:],axis=0)

    dofs_gC = 24846302 # num_nodes_gC
    nodes_gC = np.array([[13],[26],[52],[104]])
    cores_gC = nodes_gC*32
    dof_per_core_gC = dofs_gC/cores_gC

    t_avg_gC = np.array([t_gC_13_avg,t_gC_26_avg,t_gC_52_avg,t_gC_104_avg])
    t_avg_gC = np.append(nodes_gC,t_avg_gC,axis=1)
    t_avg_gC = np.append(cores_gC,t_avg_gC,axis=1)
    t_avg_gC = np.append(dof_per_core_gC,t_avg_gC,axis=1)
    t_avg_headers_gC = ['dof_per_core_gC','cores_gC','nodes_gC']
    t_avg_headers_gC = t_avg_headers_gC + th_gC_13

    linear_time_gC = t_avg_gC[0,-1]*(cores_gC[0]/cores_gC) # linear scaling

    ####################################`
    # Load g1 data
    ####################################
    file_g1_512 = root_dir+'g1.512/nrel_5mw_g1.log'
    th_g1_512,t_g1_512 = nalu.read_log(file_g1_512)
    t_g1_512_avg = np.mean(t_g1_512[-50:,:],axis=0)

    file_g1_1024 = root_dir+'g1.1024/nrel_5mw_g1.log'
    th_g1_1024,t_g1_1024 = nalu.read_log(file_g1_1024)
    t_g1_1024_avg = np.mean(t_g1_1024[-50:,:],axis=0)

    # file_g1_1536 = root_dir+'g1oarse.52/nrel_5mw_g1oarse.log'
    # th_g1_1536,t_g1_1536 = nalu.read_log(file_g1_1536)
    # t_g1_1536_avg = np.mean(t_g1_1536[500:550,:],axis=0)


    dofs_g1 = 761112205 # num_nodes_g1
    nodes_g1 = np.array([[512],[1024]])#,[1536]])
    cores_g1 = nodes_g1*32
    dof_per_core_g1 = dofs_g1/cores_g1

    t_avg_g1 = np.array([t_g1_512_avg,t_g1_1024_avg])#,t_g1_1536_avg])
    t_avg_g1 = np.append(nodes_g1,t_avg_g1,axis=1)
    t_avg_g1 = np.append(cores_g1,t_avg_g1,axis=1)
    t_avg_g1 = np.append(dof_per_core_g1,t_avg_g1,axis=1)
    t_avg_headers_g1 = ['dof_per_core_g1','cores_g1','nodes_g1']
    t_avg_headers_g1 = t_avg_headers_g1 + th_g1_512

    linear_time_g1 = t_avg_g1[0,-1]*(cores_g1[0]/cores_g1) # linear scaling

    ####################################
    ## Plots
    ####################################
    fig1 = '24.8 M Nodes (gCoarse) Timing'
    fig2 = '761.1 M Nodes (g1) Timing'
    fig3 = 'Nalu Scaling on Cori - Cores'
    fig4 = 'Nalu Scaling on Cori - DOFs per Core'

    ####################################
    # gC plotting
    ####################################
    caption_text_gC = '* NREL 5 MW on Cori Haswell noodes\n* 32 MPI ranks/node 1 OMP thread\n* Muelu solver stack with the v27.xml settings\n * 24.8 M DOF'

    plt.figure(fig1,figsize=[10,10])
    plt.title(fig1)
    for i in np.arange(1,5,1):
        plt.plot(t_gC_13[:,i],label=th_gC_13[i]+' 416 cores_gC')
        plt.plot(t_gC_26[:,i],label=th_gC_26[i]+' 832 cores_gC')
        plt.plot(t_gC_52[:,i],label=th_gC_52[i]+' 1664 cores_gC')
        plt.plot(t_gC_104[:,i],label=th_gC_104[i]+'3328 cores_gC')
    plt.legend()
    plt.xlabel('Timestep')
    plt.ylabel('Time (s)')
    plt.text(0, 100,caption_text_gC, fontsize=12)

    label = '24.8 M DOF, 32 MPI/node, 1 OMP thread, muelu v27.xml'
    plt.figure(fig3,figsize=[10,10])
    plt.title(fig3)
    plt.loglog(t_avg_gC[:,1],t_avg_gC[:,-1],'ks-',label=label)
    plt.loglog(cores_gC,linear_time_gC,'k--',label='Linear')
    plt.xlabel('Cores')
    plt.ylabel('Mean Time per Timestep (s)')
    plt.legend()

    plt.figure(fig4,figsize=[10,10])
    plt.title(fig4)
    plt.loglog(t_avg_gC[:,0],t_avg_gC[:,-1],'ks-',label=label)
    plt.loglog(dof_per_core_gC,linear_time_gC,'k--',label='linear')
    plt.xlabel('DOFs per Core')
    plt.ylabel('Mean Time per Timestep (s)')
    plt.legend()

    ####################################
    # g1 plotting
    ####################################
    caption_text_g1 = '* NREL 5 MW on Cori Haswell noodes\n* 32 MPI ranks/node 1 OMP thread\n* Muelu solver stack with the v27.xml settings\n 761.1 M DOF'

    color = 'tab:red'
    plt.figure(fig2,figsize=[10,10])
    plt.title(fig2)
    for i in np.arange(1,5,1):
        plt.plot(t_g1_512[:,i],label=th_g1_512[i]+' 16,384 cores_g1')
        plt.plot(t_g1_1024[:,i],label=th_g1_1024[i]+' 32,768 cores_g1')
        #plt.plot(t_g1_1536[:,i],label=th_g1_1536[i]+'49,152 cores_g1')

    plt.legend()
    plt.xlabel('Timestep')
    plt.ylabel('Time (s)')
    plt.text(0, 100,caption_text_g1, fontsize=12)

    label = '761.1 M DOFs, 32 MPI/node, 1 OMP thread, muelu v27.xml'
    plt.figure(fig3,figsize=[10,10])
    plt.loglog(t_avg_g1[:,1],t_avg_g1[:,-1],'s-',label=label,color=color)
    plt.loglog(cores_g1,linear_time_g1,'--',label='Linear',color=color)
    plt.xlabel('Cores')
    plt.ylabel('Mean Time per Timestep (s)')
    plt.legend()

    plt.figure(fig4,figsize=[10,10])
    plt.loglog(t_avg_g1[:,0],t_avg_g1[:,-1],'s-',label=label,color=color)
    plt.loglog(dof_per_core_g1,linear_time_g1,'--',label='linear',color=color)
    plt.xlabel('DOFs per Core')
    plt.ylabel('Mean Time per Timestep (s)')
    plt.legend()

    ####################################
    # Save plots
    ####################################
    plt.figure(fig1); plt.savefig(root_dir+fig1+'.png',dpi=400)
    plt.figure(fig2); plt.savefig(root_dir+fig2+'.png',dpi=400)
    plt.figure(fig3); plt.savefig(root_dir+fig3+'.png',dpi=400)
    plt.figure(fig4); plt.savefig(root_dir+fig4+'.png',dpi=400)
