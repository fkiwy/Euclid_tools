"""Validate tiny DR1 S3 cutouts and every actual GET range (no full mosaics).

python scripts/validate_irsa_dr1_images.py --output /tmp/irsa-dr1-images
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import astropy
import numpy as np
import s3fs
from astropy.io import fits
from astropy.wcs import WCS

from euclid_tools import irsa_tools


class AuditedS3(s3fs.S3FileSystem):
    """Observe the real anonymous S3 requests made by the production reader."""
    def __init__(self):
        super().__init__(anon=True,
                         config_kwargs={"connect_timeout": 15, "read_timeout": 60, "retries": {"max_attempts": 1}})
        self.gets = []

    async def _call_s3(self, method, *args, **kwargs):
        if method == "get_object" and not str(kwargs.get("Range", "")).startswith("bytes="):
            raise AssertionError("Refusing a full-object GET in validation")
        response = await super()._call_s3(method, *args, **kwargs)
        if method == "get_object":
            assert response["ResponseMetadata"]["HTTPStatusCode"] == 206
            self.gets.append({"key": kwargs["Key"], "range": kwargs["Range"],
                              "bytes": int(response["ContentLength"]), "content_range": response.get("ContentRange")})
        return response


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    args = cli.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    fs = AuditedS3(skip_instance_cache=True)
    cases = [("deep_VIS", "deep", "deep", "101829777", "VIS", None),
             ("deep_Y", "deep", "deep", "101829777", "Y", None),
             ("wide_survey_VIS", "wide", "survey", "102158277", "VIS", None),
             ("wide_repeat_VIS", "wide", "repeat", "102158276", "VIS",
              "s3://nasa-irsa-euclid-dr1/dr1/MER/102158276/VIS/EUC_MER_BGSUB-MOSAIC-VIS_TILE102158276-6FD3B8_20250703T031658.239336Z_00.00.fits")]
    results = []
    for name, tiling, selection, tile, band, product in cases:
        start = len(fs.gets)
        hdu = irsa_tools.retrieve_cutout(273.167,64.3687,search_radius=1,cutout_size=3,band=band,
                                        release="dr1",tiling=tiling,selection=selection,tile_id=tile,
                                        product=product,filesystem=fs)
        assert hdu is not None and hdu.data.shape == (30, 30)
        assert hdu.header["MOSCLASS"] == selection and hdu.header["TILEID"] == tile
        wcs = WCS(hdu.header)
        world = wcs.pixel_to_world(14.,15.)
        x,y = wcs.world_to_pixel(world)
        np.testing.assert_allclose([x,y],[14.,15.],atol=1e-5)
        destination = args.output / (name + ".fits")
        hdu.writeto(destination,overwrite=True)
        saved = fits.getheader(destination)
        for key in ("MAGZERO","BUNIT","PPLANID","DATASETR","PSWREL"):
            assert saved[key] == hdu.header[key]
        transfers = fs.gets[start:]
        total = sum(r["bytes"] for r in transfers)
        assert 0 < total < 1024*1024, "Tiny validation cutout unexpectedly transferred >=1 MiB"
        results.append({"name":name,"shape":list(hdu.data.shape),"source":hdu.header["SRCS3"],
                        "header":{key:hdu.header[key] for key in ("MAGZERO","BUNIT","FILTER","DATASETR","PPLANID","PIPDEFID","PSWREL","RAWHDR","CUTBYTES")},
                        "total_s3_get_bytes_including_discovery":total,"ranged_gets":transfers})
        print(name,hdu.data.shape,hdu.header["DATASETR"],total,"S3 bytes",flush=True)
    report={"utc":datetime.now(timezone.utc).isoformat(),"astropy":astropy.__version__,"s3fs":s3fs.__version__,
            "ra":273.167,"dec":64.3687,"search_radius_arcsec":1,"cutout_size_arcsec":3,"cases":results,
            "limitations":"Small VIS and Y samples only; J/H and edge/ambiguity behaviors have offline tests. Processing/calibration metadata preserved, not scientific calibration validation."}
    (args.output / "validation.json").write_text(json.dumps(report,indent=2)+'\n')


if __name__ == "__main__":
    main()
