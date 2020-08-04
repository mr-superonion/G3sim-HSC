#!/usr/bin/env python
import gc
import os
import sys
import time
import gnupg
import getpass
import configparser

import numpy as np
from astropy.table import Table, join

field_names = ['XMM',   \
              'GAMA09H',\
              'WIDE12H',\
              'GAMA15H',\
              'VVDS'   ,\
              'HECTOMAP']

gpg         =   gnupg.GPG(gnupghome="/home/xiangchong.li/.gnupg")
user        =   'xiangchong'
user_info   =   'xiangchong.li@ipmu.jp'
blinder     =   'Jim Bosch'
hsc_info    =   'hscblinder@gmail.com'

n_blind_cat =   3
dm_range    =   0.10
testing     =   True
if testing:
    user        =   'xiangchong'
    user_info   =   'xiangchong.li@ipmu.jp'
    binder  =   'xiangchong'
    hsc_info=   'xiangchong.li@ipmu.jp'

inDir   =       './.catalog_obs_reGaus_calibrated/fields/'
outDir  =       './catalog_obs_reGaus_blinded/%s/' %(user)
outDir0 =       './catalog_obs_reGaus_blinded/'
configFname =   os.path.join(outDir,'config.ini')

colname =       'i_hsmshaperegauss_derived_shear_bias_m'

def generate_dm2_values(rand_gen):
    if rand_gen<1./3.:
        m_blind=[-0.10,-0.05, 0.00]
    elif rand_gen>=1./3. and rand_gen<2./3.:
        m_blind=[-0.05, 0.00, 0.05]
    elif rand_gen>2./3.:
        m_blind=[ 0.00, 0.05, 0.10]
    return m_blind

def generate_dm1_values():
    m_blind=[]
    for ib in range(n_blind_cat):
       m_blind.append(dm_range*(2.*np.random.rand()-1))
    return m_blind

def encrypt():
    if not os.path.exists(outDir):
        os.system('mkdir -p %s' %outDir)
    # Blind random generator
    seed    =   int(time.time())
    np.random.seed(seed)

    keys    = gpg.list_keys()
    keyid1 = None
    keyid2 = None
    for k in keys:
        if user_info in k['uids'][0]:
            keyid1  =   k['fingerprint']
        if hsc_info in k['uids'][0]:
            keyid2  =   k['fingerprint']
    assert keyid1 is not None, \
            "Key fingerprint not found for %s" %(user_info)
    assert keyid2 is not None, \
            "Key fingerprint not found for %s" %(hsc_info)

    # generate dm lists
    dm1List =   generate_dm1_values()
    rand_gen=   np.random.rand()
    dm2List =   generate_dm2_values(rand_gen)

    # write the encryption info into header
    config = configparser.ConfigParser()
    config.add_section('User Blinding')
    config.set('User Blinding','keyid1', keyid1)
    config.set('User Blinding','blinder', user)
    config.add_section('HSC Blinding')
    config.set('HSC Blinding','keyid2', keyid2)
    config.set('HSC Blinding','blinder', blinder)

    seed_cry  =   gpg.encrypt('%d' %seed, keyid2).data.decode('utf8')
    config.set('HSC Blinding','seed' ,seed_cry)

    for ib in range(n_blind_cat):
        # dm1
        dm1_cry  =   gpg.encrypt('%.5f' %dm1List[ib], keyid1).data.decode('utf8')
        config.set('User Blinding','dm1-cat%d' %ib,dm1_cry)
        # dm2
        dm2_cry  =   gpg.encrypt('%.5f' %dm2List[ib], keyid2).data.decode('utf8')
        config.set('HSC Blinding','dm2-cat%d' %ib,dm2_cry)

    with open(configFname,'w') as configfile:
        config.write(configfile)

    for fieldname in field_names:
        inFname =   os.path.join(inDir,'%s_calibrated.fits' %(fieldname))
        data    =   Table.read(inFname)
        columnIn=   data[colname].copy()
        for ib in range(n_blind_cat):
            dm  =   dm2List[ib]+dm1List[ib]
            if testing and fieldname=='VVDS':
                print('dm1-cat%d: %.5f: ' %(ib,dm1List[ib]))
                print('dm2-cat%d: %.5f: ' %(ib,dm2List[ib]))
                print('dm-cat%d:  %.5f  ' %(ib,dm))
            columnOut   =   Table([columnIn+dm])
            outFname=   os.path.join(outDir,'%s_blinded_%d_m.fits' %(fieldname,ib))
            if os.path.isfile(outFname):
                print('Already have the outcome')
                continue
            columnOut.write(outFname,overwrite=True)
            del columnOut
            gc.collect()

        outNmFname=os.path.join(outDir0,'%s_no_m.fits' %(fieldname))
        if not os.path.isfile(outNmFname):
            del data[colname]
            data.write(outDir0)
        del data
        gc.collect()
    return

def decrypt():
    assert os.path.isfile(configFname),\
            'Do not have blinding info file'
    config = configparser.ConfigParser()
    config.read(configFname)
    assert config.has_section('User Blinding'),\
            "Do not have section for User Blinding"
    assert config.has_section('HSC Blinding'),\
            "Do not have section for HSC Blinding"
    keys    =   gpg.list_keys(secret=True)
    private_fingerprints = [ k['fingerprint'] for k in keys]
    keyid1  =   config.get('User Blinding','keyid1')
    if keyid1 in private_fingerprints:
        print('Decrypt User Blinding:')
        pwd = getpass.getpass(prompt="Password for %s: " %user_info)
        for ib in range(n_blind_cat):
            dm1_cry= config.get('User Blinding','dm1-cat%d' %ib)
            dec_s1 = gpg.decrypt(dm1_cry, passphrase=pwd)
            assert dec_s1.ok, "Problem during decryption: %s"%dec_s1.stderr
            dm  =   eval(dec_s1.data)
            print('dm1-cat%d: %.5f' %(ib,dm))
    else:
        print('Cannot find user\'s  fingerprint!!')

    print(private_fingerprints)
    keyid2  =   config.get('HSC Blinding','keyid2')
    if keyid2 in private_fingerprints:
        print('Decrypt HSC Blinding:')
        pwd = getpass.getpass(prompt="Password for %s: " %hsc_info)
        for ib in range(n_blind_cat):
            dm2_cry= config.get('HSC Blinding','dm2-cat%d' %ib)
            dec_s2 = gpg.decrypt(dm2_cry, passphrase=pwd)
            assert dec_s2.ok, "Problem during decryption: %s"%dec_s2.stderr
            dm  =   eval(dec_s2.data)
            print('dm2-cat%d: %.5f' %(ib,dm))
    else:
        print('Cannot find hsc\'s  fingerprint!!')
    return

def main(argv):
    if len(argv)>2:
        print('Too much arguments')
        return
    elif len(argv)<2:
        print('Arguments not enough')
        print('Please choose from (encrypt / decrypt)' )
        return
    if argv[1]=='encrypt':
        encrypt()
    elif argv[1]=='decrypt':
        decrypt()
    else:
        print('Do not support argument %s' %argv[1])
        print('Please choose from (encrypt / decrypt)' )
        return

if __name__ == "__main__":
    main(sys.argv)
