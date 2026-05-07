"""Hardware and OS information collection."""
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

from .modes import CPU_COUNT


def system_info() -> dict:
    """Collect hardware, OS, and current resource snapshot."""
    info: dict = {
        "os": f"{platform.system()} {platform.release()}",
        "arch": platform.machine(),
        "cpus": CPU_COUNT,
    }

    def _sysctl_int(key: str) -> int | None:
        r = subprocess.run(["sysctl", "-n", key], capture_output=True, text=True)
        return int(r.stdout.strip()) if r.returncode == 0 and r.stdout.strip().isdigit() else None

    def _count_cpulist(path: Path) -> int:
        n = 0
        for part in path.read_text().strip().split(","):
            if "-" in part:
                lo, hi = part.split("-", 1)
                n += int(hi) - int(lo) + 1
            elif part.strip():
                n += 1
        return n

    if sys.platform == "darwin":
        r = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                           capture_output=True, text=True)
        if r.returncode == 0:
            info["cpu"] = r.stdout.strip()
        phys = _sysctl_int("hw.physicalcpu")
        logi = _sysctl_int("hw.logicalcpu")
        if phys:
            info["physical_cores"] = phys
        if logi:
            info["logical_cores"] = logi
        nlevels = _sysctl_int("hw.nperflevels")
        if nlevels:
            cores = []
            for i in range(nlevels):
                p = _sysctl_int(f"hw.perflevel{i}.physicalcpu")
                lvl_l = _sysctl_int(f"hw.perflevel{i}.logicalcpu")
                if lvl_l:
                    label = "P" if i == 0 else "E" if i == 1 else f"L{i}"
                    entry: dict = {"type": label, "logical": lvl_l}
                    if p and p != lvl_l:
                        entry["physical"] = p
                    cores.append(entry)
            if cores:
                info["cores"] = cores
    else:
        # Linux: CPU model
        if os.path.exists("/proc/cpuinfo"):
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("model name"):
                        info["cpu"] = line.split(":", 1)[1].strip()
                        break
        cpu_dir = Path("/sys/devices/system/cpu")
        online = cpu_dir / "online"
        if online.exists():
            info["logical_cores"] = _count_cpulist(online)
        phys_ids: set[str] = set()
        for topo in cpu_dir.glob("cpu[0-9]*/topology/core_id"):
            phys_ids.add(topo.read_text().strip())
        if phys_ids:
            info["physical_cores"] = len(phys_ids)
        # Intel hybrid
        cpu_types_dir = cpu_dir / "types"
        if cpu_types_dir.is_dir():
            cores = []
            for d in sorted(cpu_types_dir.iterdir()):
                cpulist = d / "cpulist"
                if not cpulist.exists():
                    continue
                n = _count_cpulist(cpulist)
                name = d.name
                if "core" in name:
                    label = "P"
                elif "atom" in name:
                    label = "E"
                else:
                    label = name
                cores.append({"type": label, "logical": n})
            if cores:
                info["cores"] = cores
        # ARM/RISC-V
        elif (cpu_dir / "cpu0/cpu_capacity").exists():
            cap_counts: dict[int, int] = {}
            for cap_file in sorted(cpu_dir.glob("cpu[0-9]*/cpu_capacity")):
                cap = int(cap_file.read_text().strip())
                cap_counts[cap] = cap_counts.get(cap, 0) + 1
            if len(cap_counts) > 1:
                labels = ["P", "E"] + [f"L{i}" for i in range(2, 10)]
                cores = []
                for i, (cap, cnt) in enumerate(sorted(cap_counts.items(), reverse=True)):
                    cores.append({"type": labels[min(i, len(labels) - 1)], "logical": cnt})
                info["cores"] = cores

    # Total memory
    if sys.platform == "darwin":
        r = subprocess.run(["sysctl", "-n", "hw.memsize"],
                           capture_output=True, text=True)
        if r.returncode == 0:
            info["mem_total_gb"] = round(int(r.stdout.strip()) / (1 << 30), 1)
    elif os.path.exists("/proc/meminfo"):
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal"):
                    info["mem_total_gb"] = round(int(line.split()[1]) / (1 << 20), 1)
                    break

    # Available memory
    if sys.platform == "darwin":
        r = subprocess.run(["vm_stat"], capture_output=True, text=True)
        if r.returncode == 0:
            page_size = os.sysconf("SC_PAGE_SIZE") if hasattr(os, "sysconf") else 16384
            pages: dict[str, int] = {}
            for line in r.stdout.splitlines():
                if "page size of" in line:
                    page_size = int(line.split()[-2])
                elif ":" in line:
                    key, val = line.split(":", 1)
                    val = val.strip().rstrip(".")
                    if val.isdigit():
                        pages[key.strip()] = int(val)
            avail = (pages.get("Pages free", 0)
                     + pages.get("Pages inactive", 0)
                     + pages.get("Pages purgeable", 0))
            info["mem_avail_gb"] = round(avail * page_size / (1 << 30), 1)
    elif os.path.exists("/proc/meminfo"):
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable"):
                    info["mem_avail_gb"] = round(int(line.split()[1]) / (1 << 20), 1)
                    break

    # Load average
    try:
        l1, l5, l15 = os.getloadavg()
        info["load_avg"] = [round(l1, 2), round(l5, 2), round(l15, 2)]
    except OSError:
        pass

    # Memory type + Disk info
    if sys.platform == "darwin":
        try:
            r = subprocess.run(
                ["system_profiler", "SPMemoryDataType", "SPNVMeDataType", "-json"],
                capture_output=True, text=True, timeout=10)
            if r.returncode == 0:
                sp = json.loads(r.stdout)
                mem = (sp.get("SPMemoryDataType") or [{}])[0]
                mem_type = mem.get("dimm_type", "")
                mem_mfr = mem.get("dimm_manufacturer", "")
                if mem_type:
                    info["mem_type"] = f"{mem_type}" + (f" ({mem_mfr})" if mem_mfr else "")
                nvme = (sp.get("SPNVMeDataType") or [{}])[0]
                items = nvme.get("_items", [])
                if items:
                    disk = items[0]
                    info["disk"] = {
                        "model": disk.get("device_model", "?"),
                        "size": disk.get("size", "?"),
                        "protocol": "NVMe",
                        "interface": "Apple Fabric",
                    }
        except Exception:
            pass
    else:
        dmi_mem_type = Path("/sys/devices/virtual/dmi/id/memory_type")
        if dmi_mem_type.exists():
            info["mem_type"] = dmi_mem_type.read_text().strip()
        root_dev = None
        try:
            r = subprocess.run(["findmnt", "-no", "SOURCE", "/"],
                               capture_output=True, text=True, timeout=5)
            if r.returncode == 0:
                src = r.stdout.strip().split("/")[-1]
                import re
                m = re.match(r"(sd[a-z]+|nvme\d+n\d+|vd[a-z]+)", src)
                if m:
                    root_dev = m.group(1)
        except Exception:
            pass
        if root_dev:
            blk = Path(f"/sys/block/{root_dev}")
            disk_info: dict = {}
            model_f = blk / "device/model"
            if model_f.exists():
                disk_info["model"] = model_f.read_text().strip()
            rot_f = blk / "queue/rotational"
            if rot_f.exists():
                disk_info["type"] = "HDD" if rot_f.read_text().strip() == "1" else "SSD"
            if root_dev.startswith("nvme"):
                nvme_name = root_dev.split("n")[0]
                nvme_sys = Path(f"/sys/class/nvme/{nvme_name}")
                transport_f = nvme_sys / "transport"
                if transport_f.exists():
                    disk_info["protocol"] = "NVMe"
                    disk_info["interface"] = transport_f.read_text().strip().upper()
                addr_f = nvme_sys / "address"
                if addr_f.exists():
                    addr = addr_f.read_text().strip()
                    try:
                        r2 = subprocess.run(
                            ["lspci", "-vv", "-s", addr],
                            capture_output=True, text=True, timeout=5)
                        if r2.returncode == 0:
                            for line in r2.stdout.splitlines():
                                if "LnkSta:" in line and "Speed" in line:
                                    parts = line.split(",")
                                    speed = ""
                                    width = ""
                                    for p in parts:
                                        p = p.strip()
                                        if "Speed" in p:
                                            speed = p.split("Speed")[-1].strip()
                                        if "Width" in p:
                                            width = p.strip()
                                    if speed:
                                        disk_info["link_speed"] = f"{speed} {width}".strip()
                                    break
                    except Exception:
                        pass
            elif root_dev.startswith("sd"):
                disk_info["protocol"] = "SATA/SAS"
            if disk_info:
                disk_info.setdefault("size", "?")
                size_f = blk / "size"
                if size_f.exists():
                    sectors = int(size_f.read_text().strip())
                    gb = sectors * 512 / (1 << 30)
                    disk_info["size"] = f"{gb:.0f} GB"
                info["disk"] = disk_info
    return info


def format_cores(si: dict) -> str:
    """Format core topology string from system info dict."""
    phys = si.get("physical_cores")
    logi = si.get("logical_cores", si.get("cpus", "?"))
    topo = ""
    if "cores" in si:
        parts = []
        for c in si["cores"]:
            n = c.get("logical", c.get("count", "?"))
            p = c.get("physical")
            if p and p != n:
                parts.append(f"{p}{c['type']}x2")
            else:
                parts.append(f"{n}{c['type']}")
        topo = " + ".join(parts)
    if phys and phys != logi:
        if topo:
            return f"{phys}C/{logi}T ({topo})"
        return f"{phys}C/{logi}T"
    else:
        if topo:
            return f"{logi} ({topo})"
        return str(logi)
