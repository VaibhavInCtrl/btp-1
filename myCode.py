import numpy as np
import random
import matplotlib.pyplot as plt
import cv2
from datetime import datetime

# directly recycled so would be added with small delay
def Rab4Recycle(removed_particles_entering_rab4, rab4_efficiency):
    # rab4_efficiency = 0.2
    return int(np.round(removed_particles_entering_rab4*rab4_efficiency))

# sorted for recycling, so would be added with a bigger delay
def TubularRecycle(removed_particles_entering_tubular, tubular_efficiency):
    # tubular_efficiency = 0.4
    return int(np.round(removed_particles_entering_tubular*tubular_efficiency))

def Rab4Efficiency(n_part):
    if n_part < 40:
        return 0.5
    elif n_part > 100:
        return 0.2
    else:
        return -0.005 * n_part + 0.7

def RemoveRab4Efficiency(n_part):
    if n_part < 40:
        return 0.5
    elif n_part > 100:
        return 0.2
    else:
        return -0.005 * n_part + 0.7

def wildType(n_iter, mechanism_type):
    currentDate = datetime.now().date()
    currentTime = datetime.now().time()
    currentTime = str(currentTime)[:8]
    part_arr = []
    numbers = []
    rab4Number = []
    print("WILD TYPE CALLED")
    del_t = 0.01
    L=20.0
    # L is the dimension of the graph
    v = 1.0
    k = 1.0
    arr = []
    n_part = 100
    # npart is the number of molecules present on the surface
    vecx = [[] for i in range(n_part)]
    vecy = [[] for i in range(n_part)]
    # function generates the arr of molecules in the graph
    for i in range(n_part):
        x = random.random()*L
        y = random.random()*L
        vecx[i].append(x)
        vecy[i].append(y)
        arr.append((x,y))
    endo_rate = 10
    exo_rate = 10
    sigma_rate = 5
    freq = 100
    recycleFreq = 30
    recycled_particles = 0
    rab4_efficiency = 0.2

    lamRecycleFreqRab4 = 20
    lamRecycleFreqTubular = 40
    recycleFreqTubular = 40 
    recycleFreqRab4 = 20

    tubular_efficiency = 0.4
    recycleRate = 0.3
    l0=1.1 # touching distance of molecules
    img_names = []
    for t in range(n_iter):
        if t%freq==0:
            # this code runs every 100th iterations!
            # basically saves the current configuration of the surface
            fig, ax = plt.subplots()
            to_printx =[]
            to_printy =[]
            for i in range(n_part):
                to_printx.append(vecx[i][-1])
                to_printy.append(vecy[i][-1])
            # plotting the molecules location on the graph
            ax.plot(to_printx, to_printy,'o')
            fig.savefig(f"images_my/{mechanism_type}/graph_{t}.png")
            plt.close(fig)
            img_names.append(f"images_my/{mechanism_type}/graph_{t}.png")


            # REMOVING MOLECULES ENDOCYTOSIS
            nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
            # this code used to essentially make the nd (number of removed molecules less than n_part)
            while(nd>=n_part):
                print("DANGER number of removed particles are greater than particles present, bleak possibility")
                nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
            # recycled_particles = int(np.round(nd*recycleRate))
            recycleFreqTubular = np.random.poisson(lam=lamRecycleFreqTubular)
            recycleFreqRab4 = np.random.poisson(lam=lamRecycleFreqRab4)
            newValue = Rab4Efficiency(n_part)
            rab4Number.append(newValue)
            recycled_particles_rab4 = Rab4Recycle(removed_particles_entering_rab4= nd*0.5, rab4_efficiency=newValue)
            recycled_particles_tubular = TubularRecycle(removed_particles_entering_tubular=nd*0.5, tubular_efficiency=tubular_efficiency)
            # print(f"number of recycled particles \n {recycled_particles}")
            # print(f"number of recycled particles from rab4 \n {recycled_particles_rab4}")
            # print(f"number of recycled particles from tubular \n {recycled_particles_tubular}")
            
            # removed_ind has now all the indexes that are to be removed, but error here as x can be repeated
            removed_ind = set()
            while(len(removed_ind)<nd):
                x = random.randint(0, n_part-1)
                removed_ind.add(x)

            # above comment got taken care of here, but still if there was a repeat the nd would be more than the actual removed molecules
            new_indices = []
            for i in range(n_part):
                if i not in removed_ind:
                    new_indices.append(i)
            
            # replacing the actual list and positions after removing molecules
            tvecx=[]
            tvecy = []
            tarr=[]
            for ind in new_indices:
                tvecx.append(vecx[ind])
                tvecy.append(vecy[ind])
                tarr.append(arr[ind])
            vecx = tvecx
            vecy = tvecy
            arr=tarr
            # print(f"after removing molecules \n {len(arr)}")
            
            # ADDING MOLECULES EXOCYTOSIS
            na = int(np.round(random.gauss(exo_rate,sigma_rate))) - (recycled_particles_rab4 + recycled_particles_tubular)
            for i in range(na):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
            # print(f"after adding molecules \n {len(arr)}")
            # molecules are now added, so n_part is increased
            part_arr.append(n_part)

        # ADDING RECYCLES PARTICLES BACK AFTER SOME ITERATIONS I.E RECYCLED FREQ
        # elif (t+recycleFreq)%freq==0:
        #     # assuming constant endocytosis rate
        #     for i in range(recycled_particles):
        #         x = random.random()*L
        #         y = random.random()*L
        #         vecx.append([x])
        #         vecy.append([y])
        #         arr.append((x,y))
        #     n_part = len(vecx)
        #     print(f"after adding recycled molecules back \n {len(arr)}")
        elif (t+recycleFreqRab4)%freq==0:
            print(f"recycle freq for Rab 4 is \n {recycleFreqRab4}")
            # assuming constant endocytosis rate
            for i in range(recycled_particles_rab4):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
            # print(f"after adding recycled molecules from rab4 back \n {len(arr)}")
        elif (t+recycleFreqTubular)%freq==0:
            print(f"recycle freq for Tubular is \n {recycleFreqTubular}")
            # assuming constant endocytosis rate
            for i in range(recycled_particles_tubular):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
            # print(f"after adding recycled molecules from tubular back \n {len(arr)}")

        if t%10 == 0:
            numbers.append(n_part)

        # DONT CHANGE CALCULATING FORCES AND PLOTTING THE MOLECULES LOCATION EACH ITERATION
        # fxij = np.zeros(n_part)
        # fyij = np.zeros(n_part)
        # for i in range(n_part-1):
        #     for j in range(i+1,n_part):
        #         (xi,yi)=arr[i]
        #         (xj,yj)=arr[j]
        #         delxij = (xi-xj)
        #         delyij = (yi-yj)
        #         delxij = delxij - np.round(delxij/L)*L
        #         delyij = delyij - np.round(delyij/L)*L
        #         lij = np.sqrt(delxij**2.0 + delyij**2.0)
        #         if lij <= 2.0*l0:
        #             fij = -k*(lij-l0)
        #             fxij[i] += fij*delxij/abs(lij)
        #             fxij[j] -= fij*delxij/abs(lij)
        #             fyij[i] += fij*delyij/abs(lij)
        #             fyij[j] -= fij*delyij/abs(lij)
        # for i in range(n_part):
        #     (x,y)=arr[i]
        #     valid = False
        #     while(not valid):
        #         temp = random.random()
        #         theta = 2*np.pi*temp
        #         x = x + (fxij[i] + v * np.cos(theta)) * del_t
        #         y = y + (fyij[i] + v * np.sin(theta)) * del_t
        #         temp=True
        #         valid=True
        #     vecx[i].append(x)
        #     vecy[i].append(y)
        #     arr[i] = (x, y)
    total_sum = 0
    for part in range(0, len(part_arr)):
        total_sum += part_arr[part]
    print(total_sum/len(part_arr))
    plt.figure(figsize=(10,10))
    plt.plot(numbers)
    plt.ylim(bottom=0)
    plt.title('Line Graph of E-cad molecules on surface vs Time')
    plt.xlabel('Time (Every 10th iteration)')
    plt.ylabel('Number of Particles')
    plt.savefig(f"graphs/{mechanism_type}/graph_{currentDate} {currentTime}.png")
    return img_names
def Rab11(n_iter, mechanism_type):
    currentDate = datetime.now().date()
    currentTime = datetime.now().time()
    currentTime = str(currentTime)[:8]
    part_arr = []
    numbers = []
    rab4Number = []
    print("RAB11 CALLED")
    del_t = 0.01
    L=20.0
    # L is the dimension of the graph
    v = 1.0
    k = 1.0
    arr = []
    n_part = 100
    # npart is the number of molecules present on the surface
    vecx = [[] for i in range(n_part)]
    vecy = [[] for i in range(n_part)]
    # function generates the arr of molecules in the graph
    for i in range(n_part):
        x = random.random()*L
        y = random.random()*L
        vecx[i].append(x)
        vecy[i].append(y)
        arr.append((x,y))
    endo_rate = 10
    exo_rate = 10
    sigma_rate = 5
    freq = 100
    recycleFreq = 30

    lamRecycleFreqRab4 = 20
    lamRecycleFreqTubular = 40
    recycleFreqTubular = 40 
    recycleFreqRab4 = 20

    rab4_efficiency = 0.2
    tubular_efficiency = 0.4
    recycled_particles = 0
    recycleRate = 0.3
    l0=1.1 # touching distance of molecules
    img_names = []  
    for t in range(n_iter):
        if t%freq==0:
            # this code runs every 100th iterations!
            # basically saves the current configuration of the surface
            fig, ax = plt.subplots()
            to_printx =[]
            to_printy =[]
            for i in range(n_part):
                to_printx.append(vecx[i][-1])
                to_printy.append(vecy[i][-1])
            # plotting the molecules location on the graph
            ax.plot(to_printx, to_printy,'o')
            fig.savefig(f"images_my/{mechanism_type}/graph_{t}.png")
            plt.close(fig)
            img_names.append(f"images_my/{mechanism_type}/graph_{t}.png")


            # REMOVING MOLECULES ENDOCYTOSIS
            nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
            # this code used to essentially make the nd (number of removed molecules less than n_part)
            while(nd>=n_part):
                print("DANGER number of removed particles are greater than particles present, bleak possibility")
                nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
            # recycled_particles = int(np.round(nd*recycleRate))
            recycleFreqTubular = np.random.poisson(lam=lamRecycleFreqTubular)
            recycleFreqRab4 = np.random.poisson(lam=lamRecycleFreqRab4)
            newValue = Rab4Efficiency(n_part)
            rab4Number.append(newValue)
            recycled_particles_rab4 = Rab4Recycle(removed_particles_entering_rab4= nd*0.5, rab4_efficiency=newValue)
            recycled_particles_rab4_supposed = Rab4Recycle(removed_particles_entering_rab4= nd*0.5, rab4_efficiency=0.2)
            recycled_particles_tubular = TubularRecycle(removed_particles_entering_tubular=nd*0.5, tubular_efficiency=tubular_efficiency)
            recycled_particles_tubular_supposed = TubularRecycle(removed_particles_entering_tubular=nd*0.5, tubular_efficiency=0.4)
            # print(f"number of recycled particles \n {recycled_particles}")
            # print(f"number of recycled particles from rab4 \n {recycled_particles_rab4}")
            # print(f"number of recycled particles from tubular \n {recycled_particles_tubular}")
            
            # removed_ind has now all the indexes that are to be removed, but error here as x can be repeated
            removed_ind = set()
            while(len(removed_ind)<nd):
                x = random.randint(0, n_part-1)
                removed_ind.add(x)

            # above comment got taken care of here, but still if there was a repeat the nd would be more than the actual removed molecules
            new_indices = []
            for i in range(n_part):
                if i not in removed_ind:
                    new_indices.append(i)
            
            # replacing the actual list and positions after removing molecules
            tvecx=[]
            tvecy = []
            tarr=[]
            if len(new_indices) <= 50:
                rab4_efficiency = 0.3
            elif len(new_indices) >= 75:
                print("len of new_indices greater than 75")
                rab4_efficiency = 0.1
            for ind in new_indices:
                tvecx.append(vecx[ind])
                tvecy.append(vecy[ind])
                tarr.append(arr[ind])
            vecx = tvecx
            vecy = tvecy
            arr=tarr

            print(f"after removing molecules \n {len(arr)}")
            


            # ADDING MOLECULES EXOCYTOSIS
            na = int(np.round(random.gauss(exo_rate,sigma_rate))) - (recycled_particles_rab4_supposed + recycled_particles_tubular_supposed)
            for i in range(na):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
            print(f"after adding molecules \n {len(arr)}")
            # molecules are now added, so n_part is increased
            part_arr.append(n_part)

        # ADDING RECYCLES PARTICLES BACK AFTER SOME ITERATIONS I.E RECYCLED FREQ
        # elif (t+recycleFreq)%freq==0:
        #     # assuming constant endocytosis rate
        #     for i in range(recycled_particles):
        #         x = random.random()*L
        #         y = random.random()*L
        #         vecx.append([x])
        #         vecy.append([y])
        #         arr.append((x,y))
        #     n_part = len(vecx)
        #     print(f"after adding recycled molecules back \n {len(arr)}")
        elif (t+recycleFreqRab4)%freq==0:
            # assuming constant endocytosis rate
            for i in range(recycled_particles_rab4):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
            print(f"after adding recycled molecules from rab4 back \n {len(arr)}")
        # elif (t+recycleFreqTubular)%freq==0:
        #     # assuming constant endocytosis rate
        #     for i in range(recycled_particles_tubular):
        #         x = random.random()*L
        #         y = random.random()*L
        #         vecx.append([x])
        #         vecy.append([y])
        #         arr.append((x,y))
        #     n_part = len(vecx)
        #     print(f"after adding recycled molecules from tubular back \n {len(arr)}")
        if t%10 == 0:
            numbers.append(n_part)

        # DONT CHANGE CALCULATING FORCES AND PLOTTING THE MOLECULES LOCATION EACH ITERATION
        # fxij = np.zeros(n_part)
        # fyij = np.zeros(n_part)
        # for i in range(n_part-1):
        #     for j in range(i+1,n_part):
        #         (xi,yi)=arr[i]
        #         (xj,yj)=arr[j]
        #         delxij = (xi-xj)
        #         delyij = (yi-yj)
        #         delxij = delxij - np.round(delxij/L)*L
        #         delyij = delyij - np.round(delyij/L)*L
        #         lij = np.sqrt(delxij**2.0 + delyij**2.0)
        #         if lij <= 2.0*l0:
        #             fij = -k*(lij-l0)
        #             fxij[i] += fij*delxij/abs(lij)
        #             fxij[j] -= fij*delxij/abs(lij)
        #             fyij[i] += fij*delyij/abs(lij)
        #             fyij[j] -= fij*delyij/abs(lij)
        # for i in range(n_part):
        #     (x,y)=arr[i]
        #     valid = False
        #     while(not valid):
        #         temp = random.random()
        #         theta = 2*np.pi*temp
        #         x = x + (fxij[i] + v * np.cos(theta)) * del_t
        #         y = y + (fyij[i] + v * np.sin(theta)) * del_t
        #         temp=True
        #         valid=True
        #     vecx[i].append(x)
        #     vecy[i].append(y)
        #     arr[i] = (x, y)
    total_sum = 0
    for part in range(0, len(part_arr)):
        total_sum += part_arr[part]
    print(total_sum/len(part_arr))
    plt.figure(figsize=(10,10))
    plt.plot(numbers)
    plt.ylim(bottom=0)
    plt.title('Line Graph of E-cad molecules on surface vs Time')
    plt.xlabel('Time (Every 10th iteration)')
    plt.ylabel('Number of Particles')
    plt.savefig(f"graphs/{mechanism_type}/graph_{currentDate} {currentTime}.png")
    
    return img_names
def RabX1(n_iter, mechanism_type):
    currentDate = datetime.now().date()
    currentTime = datetime.now().time()
    currentTime = str(currentTime)[:8]
    part_arr = []
    numbers = []
    rab4number = []
    print("RABX1 CALLED")
    del_t = 0.01
    L=20.0
    # L is the dimension of the graph
    v = 1.0
    k = 1.0
    arr = []
    n_part = 100
    # npart is the number of molecules present on the surface
    vecx = [[] for i in range(n_part)]
    vecy = [[] for i in range(n_part)]
    # function generates the arr of molecules in the graph
    for i in range(n_part):
        x = random.random()*L
        y = random.random()*L
        vecx[i].append(x)
        vecy[i].append(y)
        arr.append((x,y))
    endo_rate = 10
    exo_rate = 10
    sigma_rate = 5
    freq = 100
    recycleFreq = 30
    lamRecycleFreqRab4 = 20
    lamRecycleFreqTubular = 40
    recycleFreqTubular = 40 
    recycleFreqRab4 = 20
    rab4_efficiency = 0.4
    tubular_efficiency = 0.08
    recycled_particles = 0
    recycleRate = 0.3
    l0=1.1 # touching distance of molecules
    img_names = []  
    for t in range(n_iter):
        if t%freq==0:
            # this code runs every 100th iterations!
            # basically saves the current configuration of the surface
            fig, ax = plt.subplots()
            to_printx =[]
            to_printy =[]
            for i in range(n_part):
                to_printx.append(vecx[i][-1])
                to_printy.append(vecy[i][-1])
            # plotting the molecules location on the graph
            ax.plot(to_printx, to_printy,'o')
            fig.savefig(f"images_my/{mechanism_type}/graph_{t}.png")
            plt.close(fig)
            img_names.append(f"images_my/{mechanism_type}/graph_{t}.png")


            # REMOVING MOLECULES ENDOCYTOSIS
            nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
            # this code used to essentially make the nd (number of removed molecules less than n_part)
            while(nd>=n_part):
                print("DANGER number of removed particles are greater than particles present, bleak possibility")
                nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
            # recycled_particles = int(np.round(nd*recycleRate))
            recycleFreqTubular = np.random.poisson(lam=lamRecycleFreqTubular)
            recycleFreqRab4 = np.random.poisson(lam=lamRecycleFreqRab4)
            newValue = RemoveRab4Efficiency(n_part)
            rab4number.append(newValue)
            recycled_particles_rab4 = Rab4Recycle(removed_particles_entering_rab4= nd*0.5, rab4_efficiency=newValue)
            recycled_particles_rab4_supposed = Rab4Recycle(removed_particles_entering_rab4= nd*0.5, rab4_efficiency=0.2)
            recycled_particles_tubular = TubularRecycle(removed_particles_entering_tubular=nd*0.5, tubular_efficiency=tubular_efficiency)
            recycled_particles_tubular_supposed = TubularRecycle(removed_particles_entering_tubular=nd*0.5, tubular_efficiency=0.4)
            # print(f"number of recycled particles \n {recycled_particles}")
            # print(f"number of recycled particles from rab4 \n {recycled_particles_rab4}")
            # print(f"number of recycled particles from tubular \n {recycled_particles_tubular}")
            
            # removed_ind has now all the indexes that are to be removed, but error here as x can be repeated
            removed_ind = set()
            while(len(removed_ind)<nd):
                x = random.randint(0, n_part-1)
                removed_ind.add(x)

            # above comment got taken care of here, but still if there was a repeat the nd would be more than the actual removed molecules
            new_indices = []
            for i in range(n_part):
                if i not in removed_ind:
                    new_indices.append(i)
            
            # replacing the actual list and positions after removing molecules
            tvecx=[]
            tvecy = []
            tarr=[]
            if len(new_indices) <= 50:
                rab4_efficiency = 0.5
            elif len(new_indices) >= 75:
                print("len of new_indices greater than 75")
                rab4_efficiency = 0.3
            for ind in new_indices:
                tvecx.append(vecx[ind])
                tvecy.append(vecy[ind])
                tarr.append(arr[ind])
            vecx = tvecx
            vecy = tvecy
            arr=tarr

            print(f"after removing molecules \n {len(arr)}")
            


            # ADDING MOLECULES EXOCYTOSIS
            na = int(np.round(random.gauss(exo_rate,sigma_rate))) - (recycled_particles_rab4_supposed + recycled_particles_tubular_supposed)
            for i in range(na):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
            print(f"after adding molecules \n {len(arr)}")
            # molecules are now added, so n_part is increased
            part_arr.append(n_part)

        # ADDING RECYCLES PARTICLES BACK AFTER SOME ITERATIONS I.E RECYCLED FREQ
        # elif (t+recycleFreq)%freq==0:
        #     # assuming constant endocytosis rate
        #     for i in range(recycled_particles):
        #         x = random.random()*L
        #         y = random.random()*L
        #         vecx.append([x])
        #         vecy.append([y])
        #         arr.append((x,y))
        #     n_part = len(vecx)
            print(f"after adding recycled molecules back \n {len(arr)}")
        elif (t+recycleFreqRab4)%freq==0:
            # print(recycleFreqRab4)
            # assuming constant endocytosis rate
            for i in range(recycled_particles_rab4):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
            print(f"after adding recycled molecules from rab4 back \n {len(arr)}")
        elif (t+recycleFreqTubular)%freq==0:
            print(recycleFreqTubular)
            # assuming constant endocytosis rate
            for i in range(recycled_particles_tubular):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
            # print(f"after adding recycled molecules from tubular back \n {len(arr)}")
        if t%10 == 0:
            numbers.append(n_part)

        # DONT CHANGE CALCULATING FORCES AND PLOTTING THE MOLECULES LOCATION EACH ITERATION
        # fxij = np.zeros(n_part)
        # fyij = np.zeros(n_part)
        # for i in range(n_part-1):
        #     for j in range(i+1,n_part):
        #         (xi,yi)=arr[i]
        #         (xj,yj)=arr[j]
        #         delxij = (xi-xj)
        #         delyij = (yi-yj)
        #         delxij = delxij - np.round(delxij/L)*L
        #         delyij = delyij - np.round(delyij/L)*L
        #         lij = np.sqrt(delxij**2.0 + delyij**2.0)
        #         if lij <= 2.0*l0:
        #             fij = -k*(lij-l0)
        #             fxij[i] += fij*delxij/abs(lij)
        #             fxij[j] -= fij*delxij/abs(lij)
        #             fyij[i] += fij*delyij/abs(lij)
        #             fyij[j] -= fij*delyij/abs(lij)
        # for i in range(n_part):
        #     (x,y)=arr[i]
        #     valid = False
        #     while(not valid):
        #         temp = random.random()
        #         theta = 2*np.pi*temp
        #         x = x + (fxij[i] + v * np.cos(theta)) * del_t
        #         y = y + (fyij[i] + v * np.sin(theta)) * del_t
        #         temp=True
        #         valid=True
        #     vecx[i].append(x)
        #     vecy[i].append(y)
        #     arr[i] = (x, y)
    total_sum = 0
    for part in range(0, len(part_arr)):
        total_sum += part_arr[part]
    print(total_sum/len(part_arr))
    plt.figure(figsize=(10,10))
    plt.plot(numbers)
    plt.ylim(bottom=0)
    plt.title('Line Graph of E-cad molecules on surface vs Time')
    plt.xlabel('Time (Every 10th iteration)')
    plt.ylabel('Number of Particles')
    plt.savefig(f"graphs/{mechanism_type}/graph_{currentDate} {currentTime}.png")
    return img_names

def modelECAD():
    currentDate = datetime.now().date()
    currentTime = datetime.now().time()
    print("Generating random value for V between 0 and 1 ")
    n_iter = int(input("Enter the number of iterations :"))
    mechanism_type = str(input("Which Type of Mechanism do you want to run? A. Wild Type B. Rab11 C. RabX1 D. All"))
    img_names = []
    
    if mechanism_type == "A" or mechanism_type == "a":
        mechanism_type = "wild_type"
        img_names = wildType(n_iter, mechanism_type)
    elif mechanism_type == "B" or mechanism_type == "b":
        mechanism_type = "rab11"
        img_names = Rab11(n_iter, mechanism_type)
    elif mechanism_type == "C" or mechanism_type == "c":
        mechanism_type = "rabx1"
        img_names = RabX1(n_iter, mechanism_type)
    elif mechanism_type == "D" or mechanism_type == "d":
        mechanism_type = "rabx1"
        rab4numberX1 = RabX1(n_iter, mechanism_type)
        mechanism_type = "rab11"
        rab4number11 = Rab11(n_iter, mechanism_type)
        mechanism_type = "wild_type"
        rab4numberWT = wildType(n_iter, mechanism_type)
        plt.figure()
        asd = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        plt.plot(rab4numberX1, label = "RabX1", linestyle = 'dashed', color = "blue")
        plt.plot(rab4number11, label = "Rab11", color = "green")
        plt.plot(rab4numberWT, label = "Wild Type", color = "grey")
        plt.ylim(bottom=0)
        plt.title('Line Graph of Rab4 Efficiency vs Time')
        plt.xlabel('Time (Every 10th iteration)')
        plt.ylabel('Rab4 Efficiency')
        plt.savefig(f"graphs/rab4/all_{currentDate} {currentTime}.png")
    print("pre-processing done!")
    codec = cv2.VideoWriter_fourcc(*"mp4v")
    out = cv2.VideoWriter(f"videos_my/{mechanism_type}/output {currentDate} {currentTime}.mp4", codec, 18, (640, 480))
    for i in range(len(img_names)):
        img = cv2.imread(img_names[i])
        out.write(img)
    plt.show()
    out.release()
    cv2.destroyAllWindows()
    print("done!")

modelECAD()