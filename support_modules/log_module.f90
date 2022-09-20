MODULE log_module
    ! The module creates the log for the material and element modules 
    ! required for debugging

    USE parameter_module, ONLY: DIRLENGTH, MSGLENGTH, NDIM_STRESS, NDIM, DP
    USE utility_module,   ONLY: PrintMatrix

    IMPLICIT NONE

    PRIVATE

    CHARACTER(len=DIRLENGTH), PUBLIC, SAVE :: logsdir

    ! Log object to report 
    TYPE, PUBLIC :: Logs
        PRIVATE
        INTEGER                    :: level   = 0
        INTEGER                    :: elem    = 0
        INTEGER                    :: inc     = 0
        CHARACTER(len=4*MSGLENGTH) :: infomsg = ''
        CHARACTER(len=MSGLENGTH)   :: warnmsg = ''
        CHARACTER(len=MSGLENGTH)   :: errmsg  = ''
    END TYPE

    PUBLIC :: set_log_level, set_log_elem, set_log_kinc, &
            &   add_log_message, write_logs 


    CONTAINS

        PURE SUBROUTINE set_log_level(log, level)
            
            TYPE(Logs), INTENT(INOUT) :: log
            INTEGER,       INTENT(IN) :: level

            log%level = level

        END SUBROUTINE set_log_level

        PURE SUBROUTINE set_log_elem(log, elem)

            TYPE(Logs), INTENT(INOUT) :: log
            INTEGER,       INTENT(IN) :: elem

            log%elem = elem

        END SUBROUTINE set_log_elem 

        PURE SUBROUTINE set_log_kinc(log, inc)
            TYPE(Logs), INTENT(INOUT) :: log
            INTEGER,       INTENT(IN) :: inc

            log%inc = inc
        END SUBROUTINE set_log_kinc

        SUBROUTINE add_log_message(log, info, warn, err)

            TYPE(Logs),              INTENT(INOUT) :: log
            CHARACTER(len=*), OPTIONAL, INTENT(IN) :: info
            CHARACTER(len=*), OPTIONAL, INTENT(IN) :: warn
            CHARACTER(len=*), OPTIONAL, INTENT(IN) :: err
            
            IF (PRESENT(info)) THEN
                IF (len(TRIM(log%infomsg)) + len(info) .GT. 4*MSGLENGTH) THEN
                    log%warnmsg = TRIM(log%warnmsg) // 'Memory warning : infomsg' // ACHAR(10)
                    CALL write_logs(log)
                    log%infomsg = TRIM(info) // ACHAR(10)
                    log%warnmsg = ''
                ELSE 
                    log%infomsg = TRIM(log%infomsg) // TRIM(info) // ACHAR(10)
                END IF
            END IF
            IF (PRESENT(warn)) log%warnmsg = TRIM(log%warnmsg) // TRIM(warn) // ACHAR(10)
            IF (PRESENT(err)) log%errmsg = TRIM(log%errmsg) // TRIM(err) // ACHAR(10)
    
        END SUBROUTINE add_log_message

        SUBROUTINE write_logs(log)
            
            USE utility_module, ONLY: newunit
            
            TYPE(Logs), INTENT(IN) :: log

            CHARACTER(len=DIRLENGTH) :: filename
            LOGICAL :: file_exist, unit_opened
            
            INTEGER :: thisunit, ierr

            thisunit = 0
            filename = ''
            
            thisunit = newunit()
            
            filename = TRIM(logsdir) // 'UEL_logs.log'
            
            INQUIRE(FILE=filename, EXIST=file_exist)

            INQUIRE(UNIT=thisunit, OPENED=unit_opened)
            
            IF(file_exist) THEN
                OPEN(thisunit, FILE=filename, STATUS='OLD', POSITION='APPEND', ACTION='WRITE')
            ELSE
                OPEN(thisunit, FILE=filename, STATUS='NEW', ACTION='WRITE', IOSTAT=ierr)
            END IF
            
            
            IF (log%level .GE. 2) THEN
                IF (LEN(TRIM(log%infomsg)) .NE. 0) THEN
                    WRITE(thisunit, '(A,I5,A,I5,A)') '[INFO, ELEM=', log%elem, ', INC=', log%inc, ']'
                    WRITE(thisunit, '(A)') log%infomsg
                END IF
                IF (LEN(TRIM(log%warnmsg)) .NE. 0) THEN
                    WRITE(thisunit, '(A,I5,A,I5,A)') '[WARNING, ELEM=', log%elem, ', INC=', log%inc, ']'
                    WRITE(thisunit, '(A)') log%warnmsg
                END IF
                IF (LEN(TRIM(log%errmsg)) .NE. 0) THEN
                    WRITE(thisunit, '(A,I5,A,I5,A)') '[ERROR, ELEM=', log%elem, ', INC=', log%inc, ']'
                    WRITE(thisunit, '(A)') log%errmsg
                END IF
            ELSE IF (log%level .EQ. 1) THEN 
                IF (LEN(TRIM(log%warnmsg)) .NE. 0) THEN
                    WRITE(thisunit, '(A,I5,A,I5,A)') '[WARNING, ELEM=', log%elem, ', INC=', log%inc, ']'
                    WRITE(thisunit, '(A)') log%warnmsg
                END IF
                IF (LEN(TRIM(log%errmsg)) .NE. 0) THEN
                    WRITE(thisunit, '(A,I5,A,I5,A)') '[ERROR, ELEM=', log%elem, ', INC=', log%inc, ']'
                    WRITE(thisunit, '(A)') log%errmsg
                END IF
            ELSE IF (log%level .EQ. 0) THEN
                IF (LEN(TRIM(log%errmsg)) .NE. 0) THEN
                    WRITE(thisunit, '(A,I5,A,I5,A)') '[ERROR, ELEM=', log%elem, ', INC=', log%inc, ']'
                    WRITE(thisunit, '(A)') log%errmsg
                END IF
            END IF

            CLOSE(thisunit)
            
        END SUBROUTINE write_logs

END MODULE log_module