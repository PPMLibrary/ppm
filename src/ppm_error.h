      !-------------------------------------------------------------------------
      !  Module       :                     ppm_error
      !-------------------------------------------------------------------------
      !
      !  Purpose      :
      !
      !  Remarks      : This is a F90 header file and not a cpp header file
      !                 since the latter would simply do replacements based on
      !                 pattern matching and does not recognize if something is
      !                 a variable name or not.
      !
      !                 Right now the default error messages are stored in a
      !                 fairly sizeable array. Maybe there is a better
      !                 solution?
      !
      !                 The default error message here should report WHAT HAS
      !                 HAPPENED. Where it happened, why and what the user can
      !                 do about it should be reported in the error message
      !                 string passed to ppm_error(.).
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_error.h,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.33  2005/06/23 18:26:19  ivos
      !  Header cosmetics.
      !
      !  Revision 1.32  2005/04/21 04:50:56  ivos
      !  Added 609.
      !
      !  Revision 1.31  2004/10/28 12:38:19  davidch
      !  Fixed numerical bug in cell lists that resulted in real
      !  particles being treated as ghosts
      !  and vice versa. The new ranking and cell list routines are
      !  supposed to be exact. All
      !  epsilons that were added to the domains in order to prevent
      !  the mentioned problems were
      !  removed since they are no longer needed.
      !  Modified Files:
      !     ppm_util_rank2d.f ppm_util_rank3d.f ppm_util_sort2d.f
      !     ppm_util_sort3d.f ppm_find_duplicates.f ppm_neighlist_clist.f
      !     ppm_error.h
      !
      !  Revision 1.30  2004/09/17 14:09:40  ivos
      !  Added 608.
      !
      !  Revision 1.29  2004/09/17 12:05:44  ivos
      !  Added 305, 704, 705.
      !
      !  Revision 1.28  2004/07/21 13:23:51  hiebers
      !  Added 805
      !
      !  Revision 1.27  2004/07/20 16:43:17  ivos
      !  Added 706.
      !
      !  Revision 1.26  2004/07/16 14:48:22  ivos
      !  Added 903.
      !
      !  Revision 1.25  2004/07/16 14:05:50  ivos
      !  Added 408 and 409.
      !
      !  Revision 1.24  2004/05/27 10:40:28  ivos
      !  added 211.
      !
      !  Revision 1.23  2004/05/26 07:37:50  ivos
      !  Added 210.
      !
      !  Revision 1.22  2004/05/13 11:39:31  ivos
      !  Added error 209.
      !
      !  Revision 1.21  2004/05/11 14:51:53  ivos
      !  Added file delete errors.
      !
      !  Revision 1.20  2004/05/06 07:28:10  ivos
      !  Added I/O error messages.
      !
      !  Revision 1.19  2004/03/31 10:46:57  ivos
      !  added 304.
      !
      !  Revision 1.18  2004/03/02 16:24:07  ivos
      !  added 703.
      !
      !  Revision 1.17  2004/02/25 13:55:27  ivos
      !  Added 607.
      !
      !  Revision 1.16  2004/02/18 12:55:34  michaebe
      !  Added 902 (routine not ready)
      !
      !  Revision 1.15  2004/02/17 16:09:27  ivos
      !  Added 605.
      !
      !  Revision 1.14  2004/02/16 12:35:26  michaebe
      !  Added entry 605: bad kick off scheme plus cosmetics.
      !
      !  Revision 1.13  2004/02/16 11:42:05  michaebe
      !  Added 901 -- 1000 and increased err_msg array size.
      !
      !  Revision 1.12  2004/02/12 14:45:40  ivos
      !  Added 304.
      !
      !  Revision 1.11  2004/02/11 14:31:44  ivos
      !  Added errors 802 through 804.
      !
      !  Revision 1.10  2004/02/10 16:56:49  hiebers
      !  added ppm_err_nofftw
      !
      !  Revision 1.9  2004/02/09 11:50:45  ivos
      !  Added 702.
      !
      !  Revision 1.8  2004/02/04 17:17:45  ivos
      !  Added 604.
      !
      !  Revision 1.7  2004/02/02 15:14:23  walther
      !  Added the error range 701 - 800.
      !
      !  Revision 1.6  2004/02/02 15:08:25  ivos
      !  Added some errors.
      !
      !  Revision 1.5  2004/01/26 12:33:02  ivos
      !  Added error 504 (used in ppm_map_part.f).
      !
      !  Revision 1.4  2004/01/23 15:18:46  ivos
      !  Added a few errors.
      !
      !  Revision 1.3  2004/01/22 13:28:49  ivos
      !  Added errors 403, 503 and 602.
      !
      !  Revision 1.2  2004/01/21 16:22:57  ivos
      !  Renamed error categories according to meeting.
      !
      !  Revision 1.1  2004/01/19 11:52:28  ivos
      !  Initial implementation.
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      CHARACTER(LEN=30), DIMENSION(1100)        :: ppm_err_mesg

      !-------------------------------------------------------------------------
      !  MPI errors:                              001 -- 100
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_err_nompi       = 001
      DATA ppm_err_mesg(001)                    /"MPI has not been initialized"/
      INTEGER, PARAMETER :: ppm_err_mpi_fail    = 002
      DATA ppm_err_mesg(002)                    /"MPI action failed"/
      INTEGER, PARAMETER :: ppm_err_mpi_term    = 003
      DATA ppm_err_mesg(003)                    /"MPI action terminated"/

      !-------------------------------------------------------------------------
      !  Memory allocation errors:                101 -- 200
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_err_alloc       = 101
      DATA ppm_err_mesg(101)                    /"Memory allocation failed"/
      INTEGER, PARAMETER :: ppm_err_dealloc     = 102
      DATA ppm_err_mesg(102)                    /"Deallocation failed"/

      !-------------------------------------------------------------------------
      !  File I/O errors:                         201 -- 300
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_err_open        = 201
      DATA ppm_err_mesg(201)                    /"Cannot open file"/
      INTEGER, PARAMETER :: ppm_err_close       = 202
      DATA ppm_err_mesg(202)                    /"Failed to close file"/
      INTEGER, PARAMETER :: ppm_err_unit_open   = 203
      DATA ppm_err_mesg(203)                    /"I/O-Unit already open"/
      INTEGER, PARAMETER :: ppm_err_no_unit     = 204
      DATA ppm_err_mesg(204)                    /"I/O-Unit not open"/
      INTEGER, PARAMETER :: ppm_err_io          = 205
      DATA ppm_err_mesg(205)                    /"I/O Error"/
      INTEGER, PARAMETER :: ppm_err_io_data     = 206
      DATA ppm_err_mesg(206)                    /"Data size missmatch"/
      INTEGER, PARAMETER :: ppm_err_file        = 207
      DATA ppm_err_mesg(207)                    /"File not found"/
      INTEGER, PARAMETER :: ppm_err_delete      = 208
      DATA ppm_err_mesg(208)                    /"Delete failed"/
      INTEGER, PARAMETER :: ppm_err_outof_units = 209
      DATA ppm_err_mesg(209)                    /"Out of I/O units"/
      INTEGER, PARAMETER :: ppm_err_data_miss   = 210
      DATA ppm_err_mesg(210)                    /"Missed data"/
      INTEGER, PARAMETER :: ppm_err_read_eor    = 211
      DATA ppm_err_mesg(211)                    /"Read past end of record"/

      !-------------------------------------------------------------------------
      !  Floating point errors:                   301 -- 400
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_err_div_zero    = 301
      DATA ppm_err_mesg(301)                    /"Division by zero"/
      INTEGER, PARAMETER :: ppm_err_sqrt_neg    = 302
      DATA ppm_err_mesg(302)                    /"Sqrt of negative number"/
      INTEGER, PARAMETER :: ppm_err_tol_warn    = 303
      DATA ppm_err_mesg(303)                    /"Unexpected tolerance"/
      INTEGER, PARAMETER :: ppm_err_bad_sum     = 304
      DATA ppm_err_mesg(304)                    /"Bad sum detected"/
      INTEGER, PARAMETER :: ppm_err_mat_singul  = 305
      DATA ppm_err_mesg(305)                    /"Matrix close to singular"/

      !-------------------------------------------------------------------------
      !  Topology and Mapping errors:             401 -- 500
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_err_nomap       = 401
      DATA ppm_err_mesg(401)                    /"No mapping defined"/
      INTEGER, PARAMETER :: ppm_err_part_range  = 402
      DATA ppm_err_mesg(402)                    /"Particles out of range found"/
      INTEGER, PARAMETER :: ppm_err_part_lost   = 403
      DATA ppm_err_mesg(403)                    /"Particles lost"/
      INTEGER, PARAMETER :: ppm_err_map_incomp  = 404
      DATA ppm_err_mesg(404)                    /"Incomplete mapping"/
      INTEGER, PARAMETER :: ppm_err_part_unass  = 405
      DATA ppm_err_mesg(405)                    /"Unassigned particles"/
      INTEGER, PARAMETER :: ppm_err_topo_missm  = 406
      DATA ppm_err_mesg(406)                    /"Particle-Topology missmatch"/
      INTEGER, PARAMETER :: ppm_err_node_number = 407
      DATA ppm_err_mesg(407)                    /"Invalid node numbering"/
      INTEGER, PARAMETER :: ppm_err_buffer_empt = 408
      DATA ppm_err_mesg(408)                    /"Buffer is empty"/
      INTEGER, PARAMETER :: ppm_err_no_topo     = 409
      DATA ppm_err_mesg(409)                    /"No such topology"/
      INTEGER, PARAMETER :: ppm_err_index_corr  = 410
      DATA ppm_err_mesg(410)                    /"Particle rank corrected"/

      !-------------------------------------------------------------------------
      !  Invalid argument:                        501 -- 600
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_err_wrong_dim   = 501
      DATA ppm_err_mesg(501)                    /"Invalid dimension"/
      INTEGER, PARAMETER :: ppm_err_wrong_prec  = 502
      DATA ppm_err_mesg(502)                    /"Invalid precision"/
      INTEGER, PARAMETER :: ppm_err_argument    = 503
      DATA ppm_err_mesg(503)                    /"Invalid argument"/
      INTEGER, PARAMETER :: ppm_err_wrong_type  = 504
      DATA ppm_err_mesg(504)                    /"Wrong data type"/

      !-------------------------------------------------------------------------
      !  Algorithmic/conceptual errors:           601 -- 700
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_err_rhs_compat  = 601
      DATA ppm_err_mesg(601)                    /"RHS violates compatibility"/
      INTEGER, PARAMETER :: ppm_err_sub_failed  = 602
      DATA ppm_err_mesg(602)                    /"Subroutine call failed"/
      INTEGER, PARAMETER :: ppm_err_bad_mesh    = 603
      DATA ppm_err_mesg(603)                    /"Bad mesh specified"/
      INTEGER, PARAMETER :: ppm_err_bad_meshop  = 604
      DATA ppm_err_mesg(604)                    /"Bad mesh operation"/
      INTEGER, PARAMETER :: ppm_err_no_data     = 605
      DATA ppm_err_mesg(605)                    /"No data"/
      INTEGER, PARAMETER :: ppm_err_bad_schko   = 606
      DATA ppm_err_mesg(606)                    /"Bad Kickoff scheme selected"/
      INTEGER, PARAMETER :: ppm_err_rev_time    = 607
      DATA ppm_err_mesg(607)                    /"Going backward in time"/
      INTEGER, PARAMETER :: ppm_err_test_fail   = 608
      DATA ppm_err_mesg(608)                    /"Test failed"/
      INTEGER, PARAMETER :: ppm_err_converge    = 609
      DATA ppm_err_mesg(609)                    /"Did not converge"/

      !-------------------------------------------------------------------------
      !  Decomposition errors:                    701 -- 800
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_err_subs_map    = 701
      DATA ppm_err_mesg(701)                    /"Sub domains not assigned"/
      INTEGER, PARAMETER :: ppm_err_mesh_miss   = 702
      DATA ppm_err_mesg(702)                    /"Mesh points missmatch"/
      INTEGER, PARAMETER :: ppm_err_subs_incomp = 703
      DATA ppm_err_mesg(703)                    /"Subs and mesh incompatible"/
      INTEGER, PARAMETER :: ppm_err_no_subs     = 704
      DATA ppm_err_mesg(704)                    /"No subdomains present"/
      INTEGER, PARAMETER :: ppm_err_few_subs    = 705
      DATA ppm_err_mesg(705)                    /"Not enough subdomains"/

      !-------------------------------------------------------------------------
      !  Library errors:                           801 -- 900
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_err_nofftw       = 801
      DATA ppm_err_mesg(801)                     /"FFTW not available"/
      INTEGER, PARAMETER :: ppm_err_nometis      = 802
      DATA ppm_err_mesg(802)                     /"METIS not available"/
      INTEGER, PARAMETER :: ppm_err_nofishpack   = 803
      DATA ppm_err_mesg(803)                     /"CRAYFISHPACK not available"/
      INTEGER, PARAMETER :: ppm_err_nohypre      = 804
      DATA ppm_err_mesg(804)                     /"HYPRE not available"/
      INTEGER, PARAMETER :: ppm_err_noMathKeisan = 805
      DATA ppm_err_mesg(805)                     /"SX MathKeisan not available"/

      !-------------------------------------------------------------------------
      !  Initialization errors:                   901 -- 1000
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_err_multipleinit = 901
      DATA ppm_err_mesg(901)                      /"Double Init attempted"/
      INTEGER, PARAMETER :: ppm_err_notready     = 902
      DATA ppm_err_mesg(902)                      /"Routine not ready"/
      INTEGER, PARAMETER :: ppm_err_ppm_noinit   = 903
      DATA ppm_err_mesg(903)                      /"PPM not initialized"/
