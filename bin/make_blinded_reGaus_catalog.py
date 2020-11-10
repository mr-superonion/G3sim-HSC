#!/usr/bin/env python
import gc
import os
import sys
import time
import gnupg
import getpass

import numpy as np
from astropy.table import Table, join

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

gpgfile     =   "/home/xiangchong.li/.gnupg"  #NOTE: The gpg file of the user
gpg         =   gnupg.GPG(gnupghome=gpgfile)

testing     =   False
if not testing:
    """
    user        =   'masamune'
    user_mail   =   'masamune.oguri@ipmu.jp'
    """
    user        =   'Wentao Luo'
    user_mail   =   'wentao.luo@ipmu.jp'
    blinder     =   'Jim Bosch'
    hsc_mail    =   'hscblinder@gmail.com'
else:
    """
    user        =   'surhud'
    user_mail   =   'surhudkicp@gmail.com'
    """
    user        =   'xiangchong'
    user_mail   =   'xiangchong.li@ipmu.jp'
    blinder      =   'xiangchong'
    hsc_mail    =   'xiangchong.li@ipmu.jp'

field_names = ['XMM',   \
              'GAMA09H',\
              'WIDE12H',\
              'GAMA15H',\
              'VVDS'   ,\
              'HECTOMAP']

def PGP_to_FITS( message):
    """
    Turns a PGP message into a string that can be stored in a FITS
    header, by removing the header and footer of the message as well
    as any \n characters
    """
    s = ""
    for st in message.split("\n"):
        if len(st) > 0 and "PGP MESSAGE" not in st and "GnuPG" not in st:
            s += st
    return s

def FITS_to_PGP(message):
    """
    Turns a string stored into an FITS comment back into a proper
    PGP message
    """
    s = "-----BEGIN PGP MESSAGE-----\n\n"
    s += message
    s += "\n-----END PGP MESSAGE-----\n"
    return s

"""
def PGP_to_FITS(message):
    #Turns a PGP message into a string that can be stored in a FITS
    #header, by removing the (\n) characters
    return message.replace('\n','=c1l=c1l=c1l=')

def FITS_to_PGP(message):
    #Turns a string stored into an FITS comment back into a proper
    #PGP message
    return message.replace('=c1l=c1l=c1l=','\n')
"""

def generate_dm2_values(rand_gen):
    if rand_gen<1./3.:
        m_blind=[-0.10,-0.05, 0.00]
    elif rand_gen>=1./3. and rand_gen<2./3.:
        m_blind=[-0.05, 0.00, 0.05]
    elif rand_gen>2./3.:
        m_blind=[ 0.00, 0.05, 0.10]
    return np.array(m_blind)

def generate_dm1_values(n_blind_cat,dm_range):
    m_blind=[]
    for ib in range(n_blind_cat):
       m_blind.append(dm_range*(2.*np.random.rand()-1))
    return m_blind

def encrypt():
    # paramters
    n_blind_cat =   1
    dm_range    =   0.10

    # directory
    inDir   =       './.catalog_obs_reGaus_calibrated/fields/'
    if testing:
        outDir  =       './catalog_obs_reGaus_blinded_m/%s-%dblind-test/' %(user.replace(' ','_'),n_blind_cat)
    else:
        outDir  =       './catalog_obs_reGaus_blinded_m/%s-%dblind/' %(user.replace(' ','_'),n_blind_cat)
    if not os.path.exists(outDir):
        os.system('mkdir -p %s' %outDir)
    outDir0 =       './catalog_obs_reGaus_no_m/'

    colname =       'i_hsmshaperegauss_derived_shear_bias_m'

    # get the key ids
    keys    = gpg.list_keys()
    keyid1 = None
    keyid2 = None
    for k in keys:
        if user_mail in k['uids'][0]:
            keyid1  =   k['fingerprint']
        if hsc_mail in k['uids'][0]:
            keyid2  =   k['fingerprint']

    # generate dm1 lists
    assert keyid1 is not None, \
            "Key fingerprint not found for %s" %(user_mail)
    dm1List =   generate_dm1_values(n_blind_cat,dm_range)
    # write the dm1 into header
    meta    =   {}
    meta.update({'ncat':n_blind_cat})
    meta.update({'user':user})
    meta.update({'uid':keyid1})

    if n_blind_cat==3:
        # generate dm2 lists
        assert keyid2 is not None, \
            "Key fingerprint not found for %s" %(hsc_mail)
        # blinding seed
        seed    =   int(time.time())
        np.random.seed(seed)
        seedtmp =   gpg.encrypt('%d' %seed, keyid2).data.decode('utf8')
        #assert '=c1l=c1l=c1l=' not in seedtmp
        seed_cry  =   PGP_to_FITS(seedtmp)
        rand_gen=   np.random.rand()
        dm2List =   generate_dm2_values(rand_gen)
        # write the dm2 into header
        meta.update({'blinder':blinder})
        meta.update({'bid':keyid2})
        meta.update({'seed':seed_cry})
        if testing: print('Random seed %d' %seed)
    else:
        dm2List =   np.zeros(n_blind_cat)


    for ib in range(n_blind_cat):
        # dm1
        dmtmp   =   gpg.encrypt('%.5f' %dm1List[ib], keyid1).data.decode('utf8')
        #assert '=c1l=c1l=c1l=' not in dmtmp
        dm1_cry =   PGP_to_FITS(dmtmp)
        meta.update({'dm1c%s' %ib:dm1_cry})
        # dm2
        if not np.all(np.abs(dm2List)<1.e-4):
            dmtmp   =  gpg.encrypt('%.5f' %dm2List[ib], keyid2).data.decode('utf8')
            #assert '=c1l=c1l=c1l=' not in dmtmp
            dm2_cry =   PGP_to_FITS(dmtmp)
            meta.update({'dm2c%s' %ib:dm2_cry})

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
            columnOut.meta.update(meta)
            outFname=   os.path.join(outDir,'%s_blinded_%d_m.fits' %(fieldname,ib))
            if os.path.isfile(outFname):
                print('Already have the blinded catalog m')
                continue
            columnOut.write(outFname,overwrite=True)
            del columnOut
            gc.collect()

        outNmFname=os.path.join(outDir0,'%s_no_m.fits' %(fieldname))
        if not os.path.isfile(outNmFname):
            del data[colname]
            data.write(outNmFname)
        del data
        gc.collect()
    return

def decrypt(header):
    n_blind_cat=    header['ncat']
    keys    =   gpg.list_keys(secret=True)
    private_fingerprints = [ k['fingerprint'] for k in keys]
    keyid1  =   header['uid']
    if keyid1 in private_fingerprints:
        print('Decrypt User Blinding:')
        pwd = getpass.getpass(prompt="Password for dm1: ")
        for ib in range(n_blind_cat):
            dm1_cry= FITS_to_PGP(header['dm1c%s' %ib])
            dec_s1 = gpg.decrypt(dm1_cry, passphrase=pwd)
            assert dec_s1.ok, "Problem during decryption: %s"%dec_s1.stderr
            dm  =   eval(dec_s1.data)
            print('dm1-cat%d: %.5f' %(ib,dm))
    else:
        print('Cannot find user\'s  fingerprint to decrypt dm1!!')

    if n_blind_cat!=3:
        print('Do not have HSC level blinding')
        return
    keyid2  =   header['bid']
    if keyid2 in private_fingerprints:
        print('\nDecrypt HSC Blinding:')
        pwd = getpass.getpass(prompt="Password for dm2: ")

        seed_cry=   FITS_to_PGP(header['seed'])
        dec_s2  =   gpg.decrypt(seed_cry, passphrase=pwd)
        assert dec_s2.ok, "Problem during decryption: %s"%dec_s2.stderr
        seed    =   eval(dec_s2.data)
        print('Random seed is: %d' %seed)
        for ib in range(n_blind_cat):
            dm2_cry= FITS_to_PGP(header['dm2c%s' %ib])
            dec_s2 = gpg.decrypt(dm2_cry, passphrase=pwd)
            assert dec_s2.ok, "Problem during decryption: %s"%dec_s2.stderr
            dm  =   eval(dec_s2.data)
            print('dm2-cat%d: %.5f' %(ib,dm))
    else:
        print('Cannot find hsc\'s  fingerprint for dm2!!')
    return

def main(argv):
    if argv[1]=='encrypt':
        if len(argv)==2:
            encrypt()
        else:
            print('Wrong number of inputs' )
    elif argv[1]=='decrypt':
        if len(argv)==3:
            if os.path.isfile(argv[2]):
                decrypt(pyfits.getheader(argv[2],1))
            else:
                print('Cannot find the fits file %s' %argv[2])
        elif len(argv)<3:
            print('Please input the relative directory of a fits file' )
        elif len(argv)>3:
            print('Wrong number of inputs' )
    else:
        print('Do not support argument %s' %argv[1])
        print('Please choose from (encrypt / decrypt)' )
        return

if __name__ == "__main__":
    main(sys.argv)
