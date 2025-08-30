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

def function_template_batch_intermediates(batch, bn, partname, declname, all_hash, all_interm_dict, method, mbpt, multiplicity, new_file, last_file, new_decl, last_decl, uses, unit_decl):

    # fix_hash = {}
    # intermediates_dict = []
    # for i in range(0, len(batch)): 
    #     fx = find_fixed_for_interm(batch[i])
    #     fix_hash[batch[i].binary_hash] = deepcopy(fx)
    #     print(batch[i], fx)
    #     minidict = {}
    #     minidict['interm'] = batch[i]
    #     minidict['int_name'] = all_hash[batch[i].binary_hash]
    #     minidict['idx_fx'] = deepcopy(fx)
    #     intermediates_dict.append(minidict)


    dct = {}
    print(method)

    if method == 'ccsd':
        name = """{partname}_interm_ccsd_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
        dname = """{declname}_interm_ccsd_pt{mbpt}""".format(declname=declname, mbpt=mbpt)
        int_app = "ss_ccsd_pt{mbpt}".format(mbpt=mbpt)
        func_name = 'ss_ccsd_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
        if multiplicity == 3:
            name = """{partname}_trip_interm_ccsd_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
            dname = """{declname}_trip_interm_ccsd_pt{mbpt}""".format(declname=declname, mbpt=mbpt)
            int_app= "tt_ccsd_pt{mbpt}".format(mbpt=mbpt)
            func_name = 'tt_ccsd_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
        if multiplicity == 13:
            name = """{partname}_so_interm_ccsd_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
            dname = """{declname}_so_interm_ccsd_pt{mbpt}""".format(declname=declname, mbpt=mbpt)
            int_app= "so_ccsd_pt_{mbpt}".format(mbpt=mbpt)
            func_name = 'so_ccsd_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
    elif method == 'cc3':
        name = """{partname}_interm_cc3_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
        dname = """{declname}_interm_cc3_pt{mbpt}""".format(declname=declname, mbpt=mbpt)
        int_app= "ss_cc3_pt{mbpt}".format(mbpt=mbpt)
        func_name = 'ss_cc3_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
        if multiplicity == 3:
            name = """{partname}_trip_interm_cc3_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
            dname = """{declname}_trip_interm_cc3_pt{mbpt}""".format(declname=declname, mbpt=mbpt)
            int_app= "tt_cc3_pt{mbpt}".format(mbpt=mbpt)
            func_name = 'tt_cc3_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
        if multiplicity == 13:
            name = """{partname}_so_interm_cc3_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
            dname = """{declname}_so_interm_cc3_pt{mbpt}""".format(declname=declname, mbpt=mbpt)
            int_app= "so_cc3_pt{mbpt}".format(mbpt=mbpt)
            func_name = 'so_cc3_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)

    dct['name'] = name
    dct['module_name'] = name
    dct['module_name_decl'] = dname

    func_name = 'pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
    dct['fname'] = func_name
    dct['mbpt'] = mbpt


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
        if '181' in minidict['int_name']:
            print(batch[i], minidict['int_name'], minidict['idx_fx'] )

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
        minidict['pt'] = mbpt
#        minidict['if'] = deepcopy(interm_fx[batch[i].binary_hash])
        intermediates_dict.append(minidict)

        print(batch[i], minidict['idx_fx'], minidict['if'], minidict['int_name'])


    name_decl = './fortran_codes/intermediates_new/decl_'+dname+'.f90'
    name_allocomp = './fortran_codes/intermediates_new/allocomp_'+name+'.f90'
    name_calls = './fortran_codes/intermediates_new/wywolania.f90'
    print(name_decl)
    print(name_allocomp)


    f_calls = open(name_calls, 'a')
    print(new_file)
    s_calls_beg = ""
    if mbpt == 0:
        if new_file == True:
            s_calls_beg = """select case(maxpt)
            """
        else:
            s_calls_beg += ""
    else:
        s_calls_beg += ""

    if new_file == True:
        s_calls_beg += """
        case({mbpt})
        """.format(**dct)
    
    s_calls = s_calls_beg + """
    call {name}_init(nocc, nactive)
    call {name}(t2, t1, s2, s1, nocc, nactive, vrdav_Rl, vrdav_Rr)
    call generate_density_exc_exc_wm(rvecl, rvecu, t2, t1, s2, s1,  &
      nocc, nactive, method, dm_wm, k1, k2, irrep0, irrep1, order, maxpt, r_idx, r_idx2,
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
    use cc3_intermediates                                                                                                                
    use arithmetic                                                                                                                  
    use s_gen                                                                                                                             
    use basis                                                                                                                      
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
    use cc3_intermediates               
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
    s_init, s_free, s = arstofort_matrix_wm(intermediates_dict, all_interm_dict, multiplicity, mbpt, int_app)

    dct['s_inner'] = s
    dct['s_init'] = s_init
    dct['s_free'] = s_free
    dct['mbpt'] = mbpt
    dct['mult'] = ""
    if multiplicity == 3:
        dct['mult'] = "triplet_"
    if multiplicity == 13:
        dct['mult'] = "so_"
    if multiplicity == 31:
        dct['mult'] = "so_left_"


    dct['multiplicity'] = multiplicity
    dct['method'] = method

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
    string = single_function_init(dct)
    string += single_function_free(dct)
    string += single_function_intermediates_wm(dct)
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

def function_template_batch(block_oo, block_ov, block_vo, block_vv, bn, partname, all_hash, method, mbpt, multiplicity, new_file, last_file, dname, s_decl, all_interm_dict):

    dct = {}
    if method == 'ccsd':
        module_name = 'ss_ccsd_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
        func_name = 'ss_ccsd_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
        int_app = "ss_ccsd_pt{mbpt}".format(mbpt=mbpt)
        if multiplicity == 3:
            module_name = 'tt_ccsd_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
            func_name = 'tt_ccsd_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
            int_app = "tt_ccsd_pt{mbpt}".format(mbpt=mbpt)
        if multiplicity == 13:
            module_name = 'so_ccsd_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
            func_name = 'so_ccsd_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
            int_app = "so_ccsd_pt{mbpt}".format(mbpt=mbpt)
#        dct['name'] = name
    elif method == 'cc3':
        module_name = 'ss_cc3_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
        int_app = "ss_cc3_pt{mbpt}".format(mbpt=mbpt)
        func_name = 'ss_cc3_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
        if multiplicity == 3:
            module_name = 'tt_cc3_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
            func_name = 'tt_cc3_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
            int_app = "tt_cc3_pt{mbpt}".format(mbpt=mbpt)
        if multiplicity == 13:
            module_name = 'so_cc3_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
            func_name = 'so_cc3_pt{mbpt}_{partname}'.format(partname=partname, mbpt=mbpt)
            int_app = "so_cc3_pt{mbpt}".format(mbpt=mbpt)


    for i in range(0, len(block_oo)):
        for j in range(0, len(block_oo[i].coefficient)):
            if 'interm' in block_oo[i].coefficient[j]:
                block_oo[i].coefficient[j] = block_oo[i].coefficient[j] + '_'+int_app

    for i in range(0, len(block_ov)):
        for j in range(0, len(block_ov[i].coefficient)):
            if 'interm' in block_ov[i].coefficient[j]:
                block_ov[i].coefficient[j] = block_ov[i].coefficient[j]+'_'+ int_app

    for i in range(0, len(block_vo)):
        for j in range(0, len(block_vo[i].coefficient)):
            if 'interm' in block_vo[i].coefficient[j]:
                block_vo[i].coefficient[j] = block_vo[i].coefficient[j]+ '_'+ int_app

    for i in range(0, len(block_vv)):
        for j in range(0, len(block_vv[i].coefficient)):
            if 'interm' in block_vv[i].coefficient[j]:
                block_vv[i].coefficient[j] = block_vv[i].coefficient[j]+ '_'+ int_app


        

    dct['name'] = func_name
    dct['module_name'] = module_name

    name_to_open = './fortran_codes/transition_excexc_new/'+module_name+'.f90'

    if new_file == True:
        f = open(name_to_open, 'w')
    else:
        f = open(name_to_open, 'a')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")
    dct['date'] = datestr

    s_beginning = """module {module_name}                                                                                                                     
    use ccsd_transformed_integrals                                                                                                    
    use cc3_intermediates                                                                                                    
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
        f.write(s_beginning)

    generate_calc_D_wm2(block_oo, 'oo', f, func_name, method, mbpt, s_decl, all_interm_dict)
    generate_calc_D_wm2(block_ov, 'ov', f, func_name, method, mbpt, s_decl, all_interm_dict)
    generate_calc_D_wm2(block_vo, 'vo', f, func_name, method, mbpt, s_decl, all_interm_dict)
    generate_calc_D_wm2(block_vv, 'vv', f, func_name, method, mbpt, s_decl, all_interm_dict)

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




