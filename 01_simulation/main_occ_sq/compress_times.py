from pycoalescence.sqlite_connection import check_sql_table_exist

def compress_times(coalescence_tree):
    if coalescence_tree.database is None:
        raise IOError("Coalescence tree has not been opened")
    species_list = [list(x) for x in coalescence_tree.get_species_list()]
    coalescence_tree._backup_species_list()
    
    for each in species_list:
        if each[6] == 1:
            each[12] = 0.0
    coalescence_tree.cursor.execute(
        "CREATE TABLE SPECIES_LIST AS SELECT * FROM SPECIES_LIST_ORIGINAL"
    )
    coalescence_tree.cursor.execute(
        "UPDATE SPECIES_LIST SET gen_added=0.0 WHERE tip == 1"
    )
    coalescence_tree.database.commit()


def double_backup(coalescence_tree):
    if coalescence_tree.database is None:
        raise IOError("Coalescence tree has not been opened")
    if check_sql_table_exist(coalescence_tree.database, "SPECIES_LIST_ORIGINAL"):
        if not check_sql_table_exist(coalescence_tree.database, "SPECIES_LIST_ORIGINAL_0"):
            coalescence_tree.cursor.execute("ALTER TABLE SPECIES_LIST_ORIGINAL RENAME TO SPECIES_LIST_ORIGINAL_0")
            coalescence_tree.database.commit()

def double_restore_backup(coalescence_tree):
    if coalescence_tree.database is None:
        raise IOError("Coalescence tree has not been opened")
    if check_sql_table_exist(coalescence_tree.database, "SPECIES_LIST_ORIGINAL"):
        coalescence_tree.revert_downsample()
    if check_sql_table_exist(coalescence_tree.database, "SPECIES_LIST_ORIGINAL_0"):
        coalescence_tree.cursor.execute("ALTER TABLE SPECIES_LIST_ORIGINAL_0 RENAME TO SPECIES_LIST_ORIGINAL")
        coalescence_tree.database.commit()
        coalescence_tree.revert_downsample()

def get_unique_times(c):
    if c.database is None:
        raise IOError("Coalescence tree has not been opened")
    times = c.cursor.execute("SELECT DISTINCT(gen_added) FROM SPECIES_LIST WHERE tip == 1").fetchall()
    return times