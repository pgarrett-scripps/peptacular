import json

from peptacular.mods.obo_parser import (
    get_entries,
    DbType,
    count_invalid_entries,
    EntryDb,
    MONOSACCHARIDES_DB,
)

# Load mono first and reset MONOSACCHARIDES_DB with the new entries
monosaccharide_entries = get_entries(
    DbType.MONOSACCHARIDES, "data/monosaccharides_updated.obo"
)

MONOSACCHARIDES_DB.reset()
MONOSACCHARIDES_DB.setup(monosaccharide_entries)


unimod_entries = get_entries(DbType.UNIMOD, "data/unimod.obo")
psi_mod_entries = get_entries(DbType.PSI_MOD, "data/psi-mod.obo")
resid_entries = get_entries(DbType.RESID, "data/psi-mod.obo")
xlmod_entries = get_entries(DbType.XLMOD, "data/xlmod.obo")
gno_entries = get_entries(DbType.GNO, "data/gno.obo")

# save entries as a json file
with open("data/entries/monosaccharide_entries.json", "w") as f:
    json.dump([e.dict() for e in monosaccharide_entries], f)

with open("data/entries/unimod_entries.json", "w") as f:
    json.dump([e.dict() for e in unimod_entries], f)

with open("data/entries/psi_mod_entries.json", "w") as f:
    json.dump([e.dict() for e in psi_mod_entries], f)

with open("data/entries/resid_entries.json", "w") as f:
    json.dump([e.dict() for e in resid_entries], f)

with open("data/entries/xlmod_entries.json", "w") as f:
    json.dump([e.dict() for e in xlmod_entries], f)

with open("data/entries/gno_entries.json", "w") as f:
    json.dump([e.dict() for e in gno_entries], f)

"""
monosaccharide_db = EntryDb(DbType.MONOSACCHARIDES, monosaccharide_entries, use_synonyms=True)
print(f'Monosaccharide None Counts (mono, ave, comp): {count_invalid_entries(monosaccharide_entries)} out of {len(monosaccharide_entries)} entries')

# add this db before updating the other dbs
monosaccharide_db.save('data/db/monosaccharide_db.pkl')


unimod_entries = get_entries(DbType.UNIMOD, 'data/unimod.obo')
psi_mod_entries = get_entries(DbType.PSI_MOD, 'data/psi-mod.obo')
resid_entries = get_entries(DbType.RESID, 'data/psi-mod.obo')
xlmod_entries = get_entries(DbType.XLMOD, 'data/xlmod.obo')
gno_entries = get_entries(DbType.GNO, 'data/gno.obo')

# print invalid entries
print(f'Unimod None Counts (mono, ave, comp): {count_invalid_entries(unimod_entries)} out of {len(unimod_entries)} entries')
print(f'PSI-MOD None Counts (mono, ave, comp): {count_invalid_entries(psi_mod_entries)} out of {len(psi_mod_entries)} entries')
print(f'RESID None Counts (mono, ave, comp): {count_invalid_entries(resid_entries)} out of {len(resid_entries)} entries')
print(f'XLMOD None Counts (mono, ave, comp): {count_invalid_entries(xlmod_entries)} out of {len(xlmod_entries)} entries')
print(f'GNO None Counts (mono, ave, comp): {count_invalid_entries(gno_entries)} out of {len(gno_entries)} entries')

# Save entries as a json file


# generate dbs
unimod_db = EntryDb(DbType.UNIMOD, unimod_entries)
psi_mod_db = EntryDb(DbType.PSI_MOD, psi_mod_entries)
resid_db = EntryDb(DbType.RESID, resid_entries)
xlmod_db = EntryDb(DbType.XLMOD, xlmod_entries)
gno_db = EntryDb(DbType.GNO, gno_entries)

# save dbs
unimod_db.save('data/db/unimod_db.pkl')
psi_mod_db.save('data/db/psi_mod_db.pkl')
resid_db.save('data/db/resid_db.pkl')
xlmod_db.save('data/db/xlmod_db.pkl')
gno_db.save('data/db/gno_db.pkl')


for db in [unimod_db, psi_mod_db, resid_db, xlmod_db, gno_db, monosaccharide_db]:
    for entry in db.id_map.values():
        print(entry.mono_mass, entry.calc_mono_mass)
"""
