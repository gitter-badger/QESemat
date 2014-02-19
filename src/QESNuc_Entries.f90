!**********************************************************************!
ENTRY QESNuc_Get_TgtName(iTarget,TgtName)
!----------------------------------------------------------------------!
    TgtName=NucName(iTarget)
    QESNuc_Get_TgtName=.true.
    return

!**********************************************************************!
ENTRY QESNuc_Get_TgtAWght(iTarget)
!----------------------------------------------------------------------!
    QESNuc_Get_TgtAWght=NucParam(2,iTarget)
    return

!**********************************************************************!
ENTRY QESNuc_Get_TgtMass(iTarget)
!----------------------------------------------------------------------!
    QESNuc_Get_TgtMass=NucMass(0,iTarget)
    return

!**********************************************************************!
ENTRY QESNuc_Get_TgtNumb(TgtName)
!----------------------------------------------------------------------!
    QESNuc_Get_TgtNumb=-999
    do Target=-1,18
        if(NucName(Target)==TgtName)then
            QESNuc_Get_TgtNumb=Target
            return                                                     !because of this |
        endif                                                          !                |
    enddo                                                              !                V
    if(QESNuc_Get_TgtNumb==-999)then                                   !we never will reach this line with QESNuc_Get_TgtNumb.ne.-999, isn't it?!
        write(*,*) 'QESNuc_Get_TgtNumb ERROR: Unknown element!'        !why not stopping just here? also -999 is needed only for this
    endif
    return

!**********************************************************************!
ENTRY QESNuc_Get_TgtNcln(iNuAnu,iTarget)                               !Caution: real function!
!----------------------------------------------------------------------!
    Z=NucParam(1,iTarget)
    A=NucParam(2,iTarget)
    if(iNuAnu==1)then
        QESNuc_Get_TgtNcln=A-Z
    else
        QESNuc_Get_TgtNcln=Z
    endif
    return

!**********************************************************************!
ENTRY QESNuc_Print_TgtParam(iTarget)
!----------------------------------------------------------------------!
    write(*,'(A40)') hf
    write(*,'(A1,1X,A9,I2,1X,A1,2X,2A2,2(A3,1X,I3,1X),A1,1X,A1)') &
    vf,'Nucleus #',iTarget,': ',adjustr(NucName(iTarget)),' (',&
    'Z =',int(NucParam(1,iTarget)),'A =',int(NucParam(2,iTarget)),')',vf
    if((iTarget.ne.0).and.(iTarget.ne.-1))then
        write(*,'(A1,11X,A27,A1)') vf,'---------------------------',vf
        write(*,'(A1,10X,"|",6X,A2,5X,"|",5X,A3,5X,A1)') vf,'Nu','Anu',vf
        write(*,'(A1,1X,A36,1X,A1)') vf,'------------------------------------',vf
        write(*,'(A1,2X,A7,2(1X,"|",1X,F11.4),1X,A1)') vf,'E_nuBin',NucParam(3,iTarget),NucParam(4,iTarget),vf
        write(*,'(A1,2X,A7,2(1X,"|",1X,F11.4),1X,A1)') vf,'P_Fermi',NucParam(5,iTarget),NucParam(6,iTarget),vf
        write(*,'(A1,2X,A7,2(1X,"|",1X,F11.4),1X,A1)') vf,'T_Fermi',NucParam(7,iTarget),NucParam(8,iTarget),vf
        write(*,'(A1,1X,A36,1X,A1)') vf,'-------------NucVolFact-------------',vf
        write(*,'(A1,2X,A7,2(1X,"|",1X,F11.8),1X,A1)') vf,'e',NucVolFact(1,1,iTarget),NucVolFact(2,1,iTarget),vf
        write(*,'(A1,2X,A7,2(1X,"|",1X,F11.8),1X,A1)') vf,'mu',NucVolFact(1,2,iTarget),NucVolFact(2,2,iTarget),vf
        write(*,'(A1,2X,A7,2(1X,"|",1X,F11.8),1X,A1)') vf,'tau',NucVolFact(1,3,iTarget),NucVolFact(2,3,iTarget),vf
    endif
    write(*,'(A40)') hf
    QESNuc_Print_TgtParam=.true.
    return
