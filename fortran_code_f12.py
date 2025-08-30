# ------------------------------------------------------                                                                                                
# Author: Aleksandra Tucholska, University of Warsaw                                                                                                     
# ------------------------------------------------------                                                                                                 
from params import *
from copy import deepcopy
import paldus_classes
from paldus_classes import *
from paldus_classes import ugg
from paldus_classes import arithmetic_string
import math
from params import virtual, occupied, general, EPSILON
from itertools import combinations
import pickle
import datetime
from eomccjac import whichvk
from eomccjac import whichvk_triplet
from eomccjac import whichvk_triplet2
from fortran_code import *




def function_template_batch_intermediates_f12(batch, bn, prefix, partname, declname, all_hash, all_interm_dict, new_file, last_file, new_decl, last_decl, uses, unit_decl):

    dct = {}

    name = """{partname}_interm_ccsd_f12""".format(partname=partname)
    dname = """{declname}_interm_ccsd_f12""".format(declname=declname)
    int_app = "{prefix}_ccsd_f12".format(prefix=prefix)
    func_name = 'ccsd_f12_{partname}'.format(partname=partname)

    dct['name'] = name
    dct['module_name'] = name
    dct['module_name_decl'] = dname

    dct['fname'] = func_name


    fix_hash = {}
    intermediates_dict = []
    for i in range(0, len(batch)):
        fx = find_fixed_for_interm(batch[i])
        fix_hash[batch[i].binary_hash] = deepcopy(fx)
        print(batch[i], fx)
        minidict = {}
        minidict['interm'] = batch[i]
        inm = all_hash[batch[i].binary_hash]
        minidict['int_name'] = all_hash[batch[i].binary_hash]+'_'+int_app
        minidict['coef_name'] = all_hash[batch[i].binary_hash]
        minidict['idx_fx'] = deepcopy(fx)
        minidict['xfx'] = all_interm_dict[inm]['xfx']
        minidict['symmetry'] = False

        # if len(minidict['xfx']) == 2:
        #     minidict['symmetry'] = True
        #     for d in minidict['xfx']:
        #         if d not in minidict['idx_fx']:
        #             print('HOLA!!!!', batch[i], minidict['idx_fx'], minidict['xfx'])
        #             sys.exit(0)
        # if minidict['symmetry'] == True:
        #     new_fx = []
        #     for d in minidict['idx_fx']:
        #         if d not in minidict['xfx']:
        #             new_fx.append(d)
        #     minidict['idx_fx'] = new_fx
        #     print('TROMBA,', fx, minidict['idx_fx'], minidict['xfx'])

                
        minidict['if'] = deepcopy(all_interm_dict[inm]['if'])
        minidict['mem'] = deepcopy(all_interm_dict[inm]['mem'])
#        minidict['if'] = deepcopy(interm_fx[batch[i].binary_hash])
        intermediates_dict.append(minidict)

        print(batch[i], minidict['idx_fx'], minidict['if'], minidict['int_name'])


    name_decl = './fortran_codes/f12-codes/decl_'+dname+'.f90'
    name_allocomp = './fortran_codes/f12-codes/allocomp_'+name+'.f90'
    name_calls = './fortran_codes/f12-codes/wywolania.f90'
    print(name_decl)
    print(name_allocomp)


    f_calls = open(name_calls, 'a')
    print(new_file)
    s_calls_beg = ""
    # if mbpt == 0:
        
    #     if new_file == True:
    #         s_calls_beg = """select case(maxpt)
    #         """
    #     else:
    #         s_calls_beg += ""
    # else:
    #     s_calls_beg += ""

    # if new_file == True:
    #     s_calls_beg += """
    #     case({mbpt})
    #     """.format(**dct)
    
    s_calls = s_calls_beg + """
    call {name}_init(nocc, nactive)
    call {name}(t2, t1, s2, s1, nocc, nactive, vrdav_Rl, vrdav_Rr)
    call generate_density_exc_exc_wm(rvecl, rvecu, t2, t1, s2, s1,  &
      nocc, nactive, dm_wm, k1, k2, irrep0, irrep1, order, maxpt, r_idx, r_idx2,
      calc_D_oo_ss_ccsd_{fname}, calc_D_oo_ss_cc3_{fname}, calc_D_oo_tt_ccsd_{fname}, calc_D_oo_tt_cc3_{fname},
      calc_D_oo_so_ccsd_{fname}, calc_D_oo_so_cc3_{fname},
      calc_D_ov_ss_ccsd_{fname}, calc_D_ov_ss_cc3_{fname}, calc_D_ov_tt_ccsd_{fname}, calc_D_ov_tt_cc3_{fname},
      calc_D_ov_so_ccsd_{fname}, calc_D_ov_so_cc3_{fname},
      calc_D_vo_ss_ccsd_{fname}, calc_D_vo_ss_cc3_{fname}, calc_D_vo_tt_ccsd_{fname}, calc_D_vo_tt_cc3_{fname},
      calc_D_vo_so_ccsd_{fname}, calc_D_vo_so_cc3_{fname},
      calc_D_vv_ss_ccsd_{fname}, calc_D_vv_ss_cc3_{fname}, calc_D_vv_tt_ccsd_{fname}, calc_D_vv_tt_cc3_{fname},
      calc_D_vv_so_ccsd_{fname}, calc_D_vv_so_cc3_{fname})
    call {name}_free()
""".format(**dct)
    
    f_calls.write(s_calls)
    f_calls.close()

    if new_decl == True:
        f_decl = open(name_decl, 'w')
    else:
        f_decl = open(name_decl, 'a')

    if new_file == True:
        f_allocomp = open(name_allocomp, 'w')
    else:
        f_allocomp = open(name_allocomp, 'a')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")
    dct['date'] = datestr

    print('DEKLARACJE')
    s_decl = ""

    for i in range(0, len(batch)):
        print(batch[i], all_hash[batch[i].binary_hash])

    for i in range(0, len(batch)):
        # if multiplicity == 1:
        #     int_name = all_hash[batch[i].binary_hash]+'_m1'
        # if multiplicity == 3:
        #     int_name = all_hash[batch[i].binary_hash]+'_m3'
        # if multiplicity == 13:
        #     int_name = all_hash[batch[i].binary_hash]+'_m13'
        int_name = all_hash[batch[i].binary_hash] + '_' + int_app 

        s_dim = ""
        for j in range(0, len(fix_hash[batch[i].binary_hash])):
            s_dim += ":, "
        s_dim = s_dim[0:len(s_dim) - 2]
        if len(fix_hash[batch[i].binary_hash]) != 0:
            s_decl += "real(F64), dimension({s_dim}), allocatable :: {int_name} \n".format(s_dim=s_dim, int_name=int_name)
        else:
            s_decl += "real(F64) :: {int_name} \n".format(s_dim=s_dim, int_name=int_name)

        
        print(batch[i], int_name)
        print(s_decl)
        print('')

    s_decl += unit_decl#"integer:: u, u1, u2, u3, u4, u5, u6, u7"

    dct['s_decl'] = s_decl
    s_decl_beginning = """module decl_{module_name_decl}
    use ccsd_transformed_integrals                                                                                                           
    use t1_transformed_int_f12
    use t1_transformed_int
    use arithmetic                                                                                                                  
    use s_gen                                                                                                                             
    use basis                                                                                                                      
    use cc_gparams
    use eom_vectors 
    implicit none
    !
    ! File generated automatically on {date}
    !
    save
    """.format(**dct)
    s_decl_decl = """
    {s_decl}
    """.format(**dct)                                                                                                                                   
    
    if new_decl == True:
        f_decl.write(s_decl_beginning)  
    f_decl.write(s_decl_decl)


    s_beginning = """module allocomp_{module_name}                                                                           
    use ccsd_transformed_integrals                     
    use t1_transformed_int_f12
    use cc_gparams
    use arithmetic                   
    use s_gen                    
    use basis                   
    use eom_vectors                   
    implicit none                   
    !                    
    ! File generated automatically on {date}                   
    !                    
    contains                    
    """.format(**dct)
    if new_file == True:
        f_allocomp.write(s_beginning)
        
    print('PRZED ARSTOFORT')
    for x in intermediates_dict:
        print(x['interm'])
    multiplicity = 1
    mbpt = 0
    s_init, s_free, s, write_decl = arstofort_intermediates_f12(intermediates_dict, all_interm_dict, int_app)

    dct['s_inner'] = s
    dct['s_init'] = s_init
    dct['s_free'] = s_free
    dct['mult'] = ""
    dct['write_decl'] = write_decl

    s_ind_lst = []
    s_ind = ""
    for i in intermediates_dict:
        for j in i['interm'].summation:
            if j not in s_ind:
                s_ind_lst.append(j)
                s_ind += "{j}, ".format(j=j)
        for j in i['idx_fx']:
            if j not in s_ind:
                s_ind_lst.append(j)
                s_ind += "{j}, ".format(j=j)

    dct['s_ind'] = s_ind[0:len(s_ind) - 2]
    dct['name_decl'] = 'decl_' + dname #'decl_'+name
    dct['uses']  = uses
    string = single_function_f12_init(dct)
    string += single_function_free(dct)
    string += single_function_intermediates_f12_SP(dct)
    f_allocomp.write(string)

    

    s_decl_end = """                                                                                                    
    end module decl_{module_name_decl}                                                                                   
    """.format(**dct)

    s_end = """                                                                                           
    end module allocomp_{module_name}                                                                              
    """.format(**dct)

    
    if last_file == True:
        f_allocomp.write(s_end)

    if last_decl == True:
        f_decl.write(s_decl_end)

    f_allocomp.close()
    f_decl.close()
    return 'decl_'+name

def function_template_batch_f12(batch_outer, prefix, bn, partname, all_hash, method, new_file, last_file, dname, s_decl, all_interm_dict):
    print('halooo', method, prefix)
    
    dct = {}
    if method == 'ccsd':
        module_name = '{prefix}_ccsd_{partname}'.format(prefix=prefix, partname=partname)
        func_name = 'calc_{prefix}_ccsd_{partname}'.format(prefix=prefix, partname=partname)
        int_app = "{prefix}_ccsd_f12".format(prefix=prefix)
        print('int apppp', int_app)
    elif method == 'cc3':
        module_name = '{prefix}_cc3_{partname}'.format(prefix=prefix, partname=partname)
        int_app = "{prefix}_cc3".format(prefix=prefix)
        func_name = '{prefix}_ccsd_{partname}'.format(prefix=prefix, partname=partname)
    print('int ap', int_app)
    for i in range(0, len(batch_outer)):
        for j in range(0, len(batch_outer[i].coefficient)):
            if 'l' in batch_outer[i].coefficient[j]:
                batch_outer[i].coefficient[j] = batch_outer[i].coefficient[j] + '_'+int_app
            print(batch_outer[i].coefficient[j], int_app)



    dct['name'] = func_name
    dct['module_name'] = module_name
    dct['decl_name'] = s_decl
    
    variables = ""
    if prefix == "t1":
        variables = "a, i"
    elif prefix == "t2":
        variables = "a, i, b, j"
    elif prefix == "tf":
        variables = "i, j, k, l"

    dct['variables'] = variables
    

    name_to_open = './fortran_codes/f12-codes/'+module_name+'.f90'

    if new_file == True:
        f = open(name_to_open, 'w')
    else:
        f = open(name_to_open, 'a')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")
    dct['date'] = datestr

    s_beginning = """module {module_name}                                                                                                                     
    use t1_transformed_int_f12
    use cc3_intermediates                                                                                                    
    use arithmetic                                                                                                   
    use s_gen                                                                                                                
    use basis     
    use cc_gparams
                                                                                                   
                                                                                                                                 
    implicit none                                                                                                                   
    !                                                                                                                                                  
    ! File generated automatically on {date}                                                                                      
    !                                                                                                               
    contains                                                                                                        
    """.format(**dct)

    if new_file == True:
        f.write(s_beginning)

    print('jjjjjestem22')
    subst_list_old = []
    subst_list_new = []
    for j in range(0, len(batch_outer)):
        subst_list_old = []
        subst_list_new = []
        print('z tego', batch_outer[j])
        for k in range(0, len(batch_outer[j].summation)):
                if ">" in batch_outer[j].summation[k]:
                    if batch_outer[j].summation[k] in virtual:
                        ii = free_idx(fortran_virtual, subst_list_new)
                        subst_list_old.append(deepcopy(batch_outer[j].summation[k]))
                        subst_list_new.append(deepcopy(ii))
                    if batch_outer[j].summation[k] in occupied:
                        ii = free_idx(fortran_occupied, subst_list_new)
                        subst_list_old.append(deepcopy(batch_outer[j].summation[k]))
                        subst_list_new.append(deepcopy(ii))
                if batch_outer[j].summation[k] in CABS:
                    ii = free_idx(fortran_CABS, subst_list_new)
                    subst_list_old.append(deepcopy(batch_outer[j].summation[k]))
                    subst_list_new.append(deepcopy(ii))
        print(subst_list_old)
        print(subst_list_new)
        for i in range(0, len(subst_list_old)):
            batch_outer[j].substitute(subst_list_old[i], subst_list_new[i])
            print("subst2", batch_outer[j])
        
    for x in batch_outer:
        print('xoxo', x)
    
    s, D, nocc, nactive = arstofort(batch_outer, 0, all_interm_dict, True, int_app)
    print(s)

    dct['s_ind'] = summation_indices(batch_outer)
    dct['D'] = D
    dct['s_inner'] = s
    string = single_function_ccsd_f12(dct)
    f.write(string)
        
    s_end = """                                                                                                                                         
    end module {module_name}                                                                                                                      
    """.format(**dct)

    if last_file == True:
        f.write(s_end)

    f.close()


def find_fixed_for_interm(e):

    fx = []
    for i in range(0, len(e.coefficient_idx)):
        for k in e.coefficient_idx[i]:
            if k not in fx and k not in e.summation:
                fx.append(k)
    return fx


def arstofort_intermediates_f12(intermediates_dict, all_interm_dict, int_app):



    interm = arithmetic_string()
    for i in range(0, len(intermediates_dict)):
        interm.append(intermediates_dict[i]['interm'])

    for i in range(0, len(interm)):
        interm[i].optimize()

    s = ""
    s_init = ""
    s_free = ""
    all_int_type = []
    write_decl = ""
    for i in range(0, len(intermediates_dict)):
        int_name = intermediates_dict[i]['int_name']
        if len(intermediates_dict[i]['if']) != 0:
            if intermediates_dict[i]['mem'] == 'ram':
                s_init += "allocate({int_name}(".format(int_name=int_name)
                for j in intermediates_dict[i]['if']:
                    if j in virt_all:
                        s_init += "nvirt0: nvirt1, "
                    elif j in occ_all:
                        s_init += "1: nocc, "
                    elif j in CABS_all:
                        s_init += "ncabs0: ncabs1, "
                s_init = s_init[0:len(s_init)-2]
                s_init += "))\n"

    for i in range(0, len(intermediates_dict)):
        if intermediates_dict[i]['mem'] == 'ram':
            int_name = intermediates_dict[i]['int_name']
            s_init += "{int_name} = zero \n".format(int_name=int_name)

    for i in range(0, len(intermediates_dict)):
        int_name = intermediates_dict[i]['int_name']
        if intermediates_dict[i]['mem'] == 'ram':
            if len(intermediates_dict[i]['if'])!= 0:
                s_free += "deallocate({int_name})\n".format(int_name=int_name)
    
    for i in range(0, len(intermediates_dict)):
        intermediates_dict[i]['interm'].optimize_f12()
                
    for i in range(0, len(intermediates_dict)):

        int_name = intermediates_dict[i]['int_name'] 
        s1 = ""
        if intermediates_dict[i]['mem'] == 'disk':
            coef_name = intermediates_dict[i]['coef_name']
            s += """open(newunit=u{coef_name}_{int_app}, file='{filename}_f12', access='stream',&
            form='unformatted', status='unknown') \n""".format(filename = all_interm_dict[coef_name]['filename'],\
                                                               coef_name=coef_name, int_app=int_app)
            write_decl += "u{coef_name}_f12, ".format(coef_name=coef_name)
            s += "allocate({coef_name}_{int_app}(".format(coef_name=coef_name, int_app=int_app)
            for j in intermediates_dict[i]['if']:
                if j in virt_all:
                    s += "nvirt0: nvirt1, "
                elif j in occ_all:
                    s += "1: nocc, "
                elif j in CABS_all:
                    s += "ncabs0: ncabs1, "                                        
            s = s[0:len(s)-2]
            s += "))\n"



        # current intermediate uses large intermediates that need to read from disk Open file and read it to temp array
        iii = 0
        for gg in range(0, len(intermediates_dict[i]['interm'].coefficient)):
            coef = intermediates_dict[i]['interm'].coefficient[gg]
            if coef in all_interm_dict:
                if all_interm_dict[coef]['mem'] == 'disk':
                    iii += 1
                    s += """open(newunit=u{coef}_{int_app}, file='{filename}_f12', access='stream',&                              
 form='unformatted', status='unknown') \n""".format(iii=iii, no = all_interm_dict[coef]['no'], \
                                                    filename = all_interm_dict[coef]['filename'], coef=coef, int_app=int_app)
                    s += "allocate({int_name}_{int_app}(".format(int_name=coef, int_app=int_app)
                    for j in all_interm_dict[coef]['if']:
                        if j in virt_all:
                            s += "nvirt0: nvirt1, "
                        elif j in occ_all:
                            s += "1: nocc, "
                        elif j in CABS_all:
                            s += "ncabs0: ncabs1, "                                        

                    s = s[0:len(s)-2]
                    s += "))\n"
                    
                    gx = ""
                    for gi in intermediates_dict[i]['interm'].coefficient_idx[gg]:
                        # if gi in virtual:
                        #     s += "do {gi} = nocc + 1, nactive \n".format(gi = gi)
                        #     nocc = True
                        #     nactive = True
                        # if gi in occupied:
                        #     s += "do {gi} = 1, nocc \n".format(gi = gi)
                        #     nocc = True
                        gx += "{gi}, ".format(gi=gi)
                    gx = gx[0:len(gx)-2]
                    
                    s += 'read (u{coef}_{int_app}) {coef}_{int_app}\n'.format(coef = coef, gx=gx, int_app=int_app)
                    # for gi in range(0, len(intermediates_dict[i]['interm'].coefficient_idx[gg])):
                    #     s += "end do \n"            

        subst_list_old = []
        subst_list_new = []
                    
        for j in intermediates_dict[i]['interm'].summation:
            print('sprawdzam indeks', j)
            if ">" in j:
                if j in virtual:
                    ii = free_idx(fortran_virtual, subst_list_new)
                    subst_list_old.append(deepcopy(j))
                    subst_list_new.append(deepcopy(ii))
                    idx = intermediates_dict[i]['interm'].summation.index(j)
                    intermediates_dict[i]['interm'].summation[idx] = ii
                if j in occupied:
                    print('TAK IN OCCC')
                    ii = free_idx(fortran_occupied, subst_list_new)
                    subst_list_old.append(deepcopy(j))
                    subst_list_new.append(deepcopy(ii))
                    idx = intermediates_dict[i]['interm'].summation.index(j)
                    intermediates_dict[i]['interm'].summation[idx] = ii
            if j in CABS:
                ii = free_idx(fortran_CABS, subst_list_new)
                subst_list_old.append(deepcopy(j))
                subst_list_new.append(deepcopy(ii))
                idx = intermediates_dict[i]['interm'].summation.index(j)
                intermediates_dict[i]['interm'].summation[idx] = ii

        for j in intermediates_dict[i]['idx_fx']:
            print('sprawdzam indeks', j)
            if ">" in j:
                if j in virtual:
                    ii = free_idx(fortran_virtual, subst_list_new)
                    subst_list_old.append(deepcopy(j))
                    subst_list_new.append(deepcopy(ii))
                    idx = intermediates_dict[i]['idx_fx'].index(j)
                    intermediates_dict[i]['idx_fx'][idx] = ii
                if j in occupied:
                    print('TAK IN OCCC')
                    ii = free_idx(fortran_occupied, subst_list_new)
                    subst_list_old.append(deepcopy(j))
                    subst_list_new.append(deepcopy(ii))
                    idx = intermediates_dict[i]['idx_fx'].index(j)
                    intermediates_dict[i]['idx_fx'][idx] = ii
            if j in CABS:
                ii = free_idx(fortran_CABS, subst_list_new)
                subst_list_old.append(deepcopy(j))
                subst_list_new.append(deepcopy(ii))
                idx = intermediates_dict[i]['idx_fx'].index(j)
                intermediates_dict[i]['idx_fx'][idx] = ii


        print("0subst1", intermediates_dict[i]['interm'], interm[i])
        print(subst_list_old)
        print(subst_list_new)
        for j in range(0, len(subst_list_old)):
            print('jjj', j, intermediates_dict[i]['if'])
            interm[i].substitute(subst_list_old[j], subst_list_new[j])
            if subst_list_old[j] in intermediates_dict[i]['if']:
                idx = intermediates_dict[i]['if'].index(subst_list_old[j])
                intermediates_dict[i]['if'][idx] = subst_list_new[j]
        print("0subst2", intermediates_dict[i]['interm'], interm[i])

        for j in intermediates_dict[i]['interm'].summation:
            s1 += "{j}, ".format(j=j)
        for j in intermediates_dict[i]['idx_fx']:
            s1 += "{j}, ".format(j=j)
        s1 += "sum"
        if (len(intermediates_dict[i]['idx_fx'])!= 0):
            s += "!$omp parallel private({s1})& \n".format(s1=s1)
            s += "!$omp default(shared)\n"
            s += "!$omp do collapse({n})\n".format(n=len(intermediates_dict[i]['idx_fx']))
        for j in intermediates_dict[i]['idx_fx']:
            if j in virt_all:
                s += "do {j} = nvirt0, nvirt1 \n".format(j = j)
                nocc = True
                nactive = True
            if j in occ_all:
                s += "do {j} = 1, nocc \n".format(j = j)
                nocc = True
            if j in CABS_all:
                s += "do {j} = ncabs0, ncabs1 \n".format(j = j)
                ncabs = True
                nocc = True



        s += "sum = zero \n"
        for j in interm[i].summation:
            if j in virt_all:
                s += "do {j} = nocc + 1, nactive \n".format(j = j)
                nocc = True
                nactive = True
            if j in occ_all:
                s += "do {j} = 1, nocc \n".format(j = j)
                nocc = True
            if j in CABS_all:
                s += "do {j} = ncabs0, ncabs1 \n".format(j = j)
                ncabs = True
                nocc = True


        s += "sum = sum + "  
        for j in range(0, len(interm[i].coefficient)):

            if interm[i].coefficient[j] == CC_AMPLITUDE:
                ln, at = interm[i].amplitude_type(j)
                if ln == '3' :
                    temp = "(nocc, nactive, "
                else:
                    temp = "("
                for k in interm[i].coefficient_idx[j]:
                    temp = temp + str(k) + ","
                temp = temp[0:len(temp)-1]
                if (j + 1) < len(interm[i].coefficient):
                    s += "{at}{temp}) * ".format(at = at, temp = temp)
                else:
                    s += "{at}{temp})".format(at = at, temp = temp)


            elif interm[i].coefficient[j] == S_AMPLITUDE:
                ln, at = interm[i].amplitude_type(j)
                if ln == '3' :
                    temp = "(nocc, nactive, "
                else:
                    temp = "("
                for k in interm[i].coefficient_idx[j]:
                    temp = temp + str(k) + ","
                temp = temp[0:len(temp)-1]
                if (j + 1) < len(interm[i].coefficient):
                    s += "{at}{temp}) * ".format(at = at, temp = temp)
                else:
                    s += "{at}{temp})".format(at = at, temp = temp)

            elif interm[i].coefficient[j] == INTERM_Ft_F12:
                at = "Ft_interm"
                temp = "("
                for k in interm[i].coefficient_idx[j]:
                    temp = temp + str(k) + ","
                temp = temp[0:len(temp)-1]
                if (j + 1) < len(interm[i].coefficient):
                    s += "{at}{temp}) * ".format(at = at, temp = temp)
                else:
                    s += "{at}{temp})".format(at = at, temp = temp)

            elif interm[i].coefficient[j] == INTERM_V_F12:
                g_type = interm[i].int_type("vvf", j)
                print('gr1', g_type)
                if(g_type not in all_int_type):
                    all_int_type.append(g_type)

                idx = interm[i].coefficient_idx[j]
                if (j + 1) < len(interm[i].coefficient):
                    s += ("{g_type}({i0}, {i1}, {i2}, {i3}) * ".format(g_type = g_type, 
                                                                       i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                else :
                    s += ("{g_type}({i0}, {i1}, {i2}, {i3})".format(g_type = g_type, 
                                                                    i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
            elif interm[i].coefficient[j] == INTERM_X_F12:

                g_type = interm[i].int_type("xff", j)
                print('gr2', g_type)
                if(g_type not in all_int_type):
                    all_int_type.append(g_type)

                idx = interm[i].coefficient_idx[j]
                if (j + 1) < len(interm[i].coefficient):
                    s += ("{g_type}({i0}, {i1}, {i2}, {i3}) * ".format(g_type = g_type, 
                                                    i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                else :
                    s += ("{g_type}({i0}, {i1}, {i2}, {i3})".format(g_type = g_type, 
                                                i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))

                    
            elif interm[i].coefficient[j] == INTERM_B_F12:

                g_type = interm[i].int_type("bff", j)
                print('gr3', g_type)
                if(g_type not in all_int_type):
                    all_int_type.append(g_type)

                idx = interm[i].coefficient_idx[j]
                if (j + 1) < len(interm[i].coefficient):
                    s += ("{g_type}({i0}, {i1}, {i2}, {i3}) * ".format(g_type = g_type, 
                                                    i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                else :
                    s += ("{g_type}({i0}, {i1}, {i2}, {i3})".format(g_type = g_type, 
                                                i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                
                
            elif interm[i].coefficient[j] == INTERM_P_F12:
                g_type = interm[i].int_type("pfvf", j)
                print('gr4', g_type)
                if(g_type not in all_int_type):
                    all_int_type.append(g_type)

                idx = interm[i].coefficient_idx[j]
                if (j + 1) < len(interm[i].coefficient):
                    s += ("{g_type}({i0}, {i1}, {i2}, {i3}) * ".format(g_type = g_type, 
                                                    i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                else :
                    s += ("{g_type}({i0}, {i1}, {i2}, {i3})".format(g_type = g_type, 
                                                i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))

            elif interm[i].coefficient[j] == F12_AMPLITUDE:
                at = "t2f"
                temp = "("
                for k in interm[i].coefficient_idx[j]:
                    temp = temp + str(k) + ","
                temp = temp[0:len(temp)-1]
                if (j + 1) < len(interm[i].coefficient):
                    s += "{at}{temp}) * ".format(at = at, temp = temp)
                else:
                    s += "{at}{temp})".format(at = at, temp = temp)

            elif interm[i].coefficient[j] == TWOEL_INT:
                g_type = interm[i].int_type("", j)
                print('gr5', g_type)
                if(g_type not in all_int_type):
                    all_int_type.append(g_type)

                idx = interm[i].coefficient_idx[j]
    
                if (j + 1) < len(interm[i].coefficient):
                    s += ("{g_type}({i0}, {i1}, {i2}, {i3}) * ".format(g_type = g_type, 
                                                                       i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                else :
                    s += ("{g_type}({i0}, {i1}, {i2}, {i3})".format(g_type = g_type, 
                                                                    i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))

            elif interm[i].coefficient[j] == F12_TWOEL:
                print('takjestemo')
                g_type = interm[i].int_type("ff", j)
                print('gr6', g_type)
                if(g_type not in all_int_type):
                    all_int_type.append(g_type)

                idx = interm[i].coefficient_idx[j]
    
                if (j + 1) < len(interm[i].coefficient):
                    s += ("{g_type}({i0}, {i1}, {i2}, {i3}) * ".format(g_type = g_type, 
                                                                       i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                else :
                    s += ("{g_type}({i0}, {i1}, {i2}, {i3})".format(g_type = g_type, 
                                                                    i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))


            elif interm[i].coefficient[j] == TWOEL_INT_TRANS:
                
                v_type = interm[i].int_type("t", j)
                print('gr7', v_type)
                if(v_type not in all_int_type):
                    all_int_type.append(v_type)

                idx = interm[i].coefficient_idx[j]

                if (j + 1) < len(interm[i].coefficient):
                    s += ("{v_type}({i0}, {i1}, {i2}, {i3}) * ".format(v_type = v_type, 
                                                                i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                else :

                    s += ("{v_type}({i0}, {i1}, {i2}, {i3})".format(v_type = v_type, 
                                                                    i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))

            elif interm[i].coefficient[j] == BARENUCL_HAM_TRANS:
                idx = interm[i].coefficient_idx[j]
                k_type = interm[i].int1_type_trans(j)

                if (j + 1) < len(interm[i].coefficient):
                    s += "{k_type}({i0}, {i1}) * ".format(k_type = k_type, i0 = idx[0], i1 = idx[1])
                else :
                    s += "{k_type}({i0}, {i1})".format(k_type = k_type, i0 = idx[0], i1 = idx[1])


            elif interm[i].coefficient[j] == FOCK_MATRIX or interm[i].coefficient[j] == FOCK_MATRIX_TRANS:
                idx = interm[i].coefficient_idx[j]
                if idx[0] == idx[1]:
                    if (j + 1) < len(interm[i].coefficient):
                        s += "eorb({i0}) * ".format(i0 = idx[0])
                    else :
                        s += "eorb({i0})".format(i0 = idx[0])


            else:
                temp = "("
                for k in interm[i].coefficient_idx[j]:
                    temp = temp + str(k) + ","
                temp = temp[0:len(temp)-1]
                if (j + 1) < len(interm[i].coefficient):
                    s += "{coef}_{int_app}{temp}) * ".format(coef = interm[i].coefficient[j], int_app = int_app, temp = temp)
                else:
                    s += "{coef}_{int_app}{temp})".format(coef = interm[i].coefficient[j], int_app = int_app, temp = temp)


                
        s += "\n"

        for j in range(0, len(interm[i].summation)):
            s += "end do \n"

        if len(intermediates_dict[i]['if']) != 0:
            print(1, int_name)

            m = "{int_name}(".format(int_name=int_name)
            for k in intermediates_dict[i]['if']:
                m += "{k}, ".format(k=k)
            m = m[0:len(m)-2]
            m += ")"
        else:
            print(2, int_name)
            m = "{int_name}".format(int_name=int_name)


        m += " = + sum \n"
        s += m
       

        for j in range(0, len(intermediates_dict[i]['if'])):
            s += "end do \n"
        if (len(intermediates_dict[i]['if'])!= 0):
            s += "!$omp end do nowait \n"
            s += "!$omp end parallel \n"
            s += "\n"

        if intermediates_dict[i]['mem'] == 'disk':


            coef_name = intermediates_dict[i]['coef_name']
            gx = ""
            for rr in intermediates_dict[i]['if']:
                gx += rr + ", "
            gx = gx[0:len(gx)-2]
            s+= "write(u{coef_name}_{int_app}) {coef_name}_{int_app} \n".format(coef_name=coef_name, gx=gx, int_app = int_app)
            s += "deallocate({coef_name}_{int_app}) \n".format(coef_name=coef_name, int_app=int_app)
            s += 'close(u{coef_name}_{int_app}) \n'.format(coef_name=coef_name, int_app=int_app)

        iii = 0
        for gg in range(0, len(intermediates_dict[i]['interm'].coefficient)):
            coef = intermediates_dict[i]['interm'].coefficient[gg]
            if coef in all_interm_dict:
                if all_interm_dict[coef]['mem'] == 'disk':
                    iii += 1
                    s += """close(u{coef}_{int_app}) \n""".format(coef=coef, int_app=int_app)
                    s += "deallocate({coef}_{int_app}) \n".format(coef=coef, int_app=int_app)

    print('all int types')
    for x in all_int_type:
        print(x)
    print('')

    # if write_decl[len(write_decl)-1]==",":
    #     write_decl = write_decl[0:len(write_decl)-2]

    return s_init, s_free, s, write_decl
