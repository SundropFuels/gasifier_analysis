USE gas_unit_test;

DROP TABLE IF EXISTS dataframe_pd_test_table;

CREATE TABLE dataframe_pd_test_table
(
    A float,
    B float,
    C float,
    cheetah float
);

INSERT INTO dataframe_pd_test_table VALUES(1.2,4.6,3.6,2.6);
INSERT INTO dataframe_pd_test_table VALUES(3.1,7.0,8.0,9.2);
INSERT INTO dataframe_pd_test_table VALUES(1.1,7.3,2.5,1.1);


DROP TABLE IF EXISTS dataframe_upload_test_table;

CREATE TABLE dataframe_upload_test_table
(
    A float,
    B float,
    C float,
    cheetah float
);






