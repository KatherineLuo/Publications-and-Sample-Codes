SELECT spc.specimen_id, spc.barcode, spc.req_specimen_id, spc.sp_ab_id, spc.specimen_class, spc.specimen_type, spc.parent_specimen_id, 
       spc.specimen_characteristics_id, spc.char_id, spc.tissue_side, spc.tissue_site, spc.collection_protocol_id, spc.scg_id, spc.scg_name, spc.collection_point_label, 
       spc.collection_protocol_reg_id, spc.coll_prot_reg_id, spc.protocol_participant_id, spc.registration_date, spc.participant_id,  
       pos.ab_pos_id, pos.position_dimension_one, pos.position_dimension_two, pos.sp_pos_id, pos.specimen_id, pos.container_id
FROM (
  SELECT sp.specimen_id, sp.barcode, sp.specimen_collection_group_id, sp.req_specimen_id, 
         sp.sp_ab_id, sp.specimen_class, sp.specimen_type, sp.parent_specimen_id, sp.specimen_characteristics_id, sp.char_id, sp.tissue_side, sp.tissue_site,
         f.collection_protocol_id, f.scg_id, f.scg_name, f.collection_point_label, f.collection_protocol_reg_id, f.coll_prot_reg_id, f.protocol_participant_id, f.registration_date, f.participant_id
  FROM (SELECT b.collection_protocol_id, A.scg_id, A.NAME scg_name, A.collection_point_label, A.collection_protocol_reg_id, 
               b.coll_prot_reg_id, b.protocol_participant_id, b.registration_date, b.participant_id  
        FROM 
           (SELECT a1.scg_id, a1.NAME, a1.collection_protocol_reg_id, b1.collection_point_label
            FROM ( SELECT IDENTIFIER scg_id, NAME, collection_protocol_reg_id, collection_protocol_event_id
                   FROM catissue.catissue_specimen_coll_group 
                   WHERE collection_protocol_event_id IN (SELECT IDENTIFIER 
                                                          FROM catissue.catissue_coll_prot_event
                                                          WHERE Collection_protocol_id = 115)
                  ) a1,
           
                   CATISSUE.catissue_coll_prot_event b1
           
            WHERE a1.collection_protocol_event_id = b1.IDENTIFIER) A,
 
           (SELECT IDENTIFIER coll_prot_reg_id, protocol_participant_id, registration_date, participant_id, collection_protocol_id 
            FROM CATISSUE.catissue_coll_prot_reg 
            WHERE collection_protocol_id = 115) b

          WHERE A.collection_protocol_reg_id(+) = b.coll_prot_reg_id ) f,
        
       (SELECT y.IDENTIFIER specimen_id, y.barcode, y.specimen_collection_group_id, y.req_specimen_id, 
               w.IDENTIFIER sp_ab_id, w.specimen_class, w.specimen_type, w.parent_specimen_id, w.specimen_characteristics_id, z.identifier char_id, z.tissue_side, z.tissue_site
       FROM catissue.caTissue_specimen y, catissue.caTissue_abstract_specimen w, caTissue.catissue_specimen_char z
       WHERE y.identifier(+) = w.identifier AND w.specimen_characteristics_id = z.identifier ) sp
        
  WHERE f.scg_id = sp.specimen_collection_group_id(+)  ) spc,
  
  (SELECT ab_pos.IDENTIFIER ab_pos_id, ab_pos.position_dimension_one, ab_pos.position_dimension_two, sp_pos.IDENTIFIER sp_pos_id, sp_pos.specimen_id, sp_pos.container_id
   FROM caTissue.caTissue_abstract_position ab_pos,
        caTissue.caTissue_specimen_position sp_pos
   WHERE ab_pos.IDENTIFIER = sp_pos.IDENTIFIER) pos
   
WHERE spc.specimen_id = pos.specimen_id(+) ;

