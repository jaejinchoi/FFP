### [Note] 2v.3.1
* From version 2v.3.0, only binary FFP will be continued and supported (for text-based FFP the latest version is 2v.2.1)
* FFP_bin's value presentation is Long Long (for feature count) or Double (for feature frequency).  
  A variable byte definition dependent on platforms because it is byte based.    
* JSD_matrix_bin's value presentation is a point below 8 decimal places (%.8g)  
* Option [-v] shows a program profile  
* 2v.3.0 fully supports user-defined option [-k] 


2020-2-24  
* In JSD_matrix calculation, a valid-time of 'q_f_buf' may cause adverse consequences and so there was a code rearrangement.  

2020-1-13  
* Unadvised to input large genome files. It shown to loss/malfunction during zlib compression.  

