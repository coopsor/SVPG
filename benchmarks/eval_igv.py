import socket
import time
import os
from PIL import Image, ImageOps
import collections


def connect_igv(host='127.0.0.1', port=60151):
    """Connect to IGV."""
    try:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect((host, port))
        print("Successfully connected to IGV")
        return sock
    except ConnectionRefusedError:
        print("Failed to connect to IGV. Make sure IGV is running and listening on port 60151.")
        raise
    except Exception as e:
        print(f"Error occurred while connecting to IGV: {str(e)}")
        raise


def send_command(sock, command):
    """Send command to IGV and wait for response."""
    try:
        print(f"Sending command: {command}")
        sock.send(command.encode() + b'\n')

        # Wait time depending on the command type
        if 'load' in command:
            print("Waiting for file to load... (5 seconds)")
            time.sleep(5)
        elif 'genome' in command:
            print("Waiting for genome to load... (2 seconds)")
            time.sleep(2)
        elif 'goto' in command:
            print("Waiting to jump to location... (1 second)")
            time.sleep(1)
        elif 'snapshot' in command:
            print("Waiting for snapshot... (1 second)")
            time.sleep(1)
        else:
            time.sleep(1)

        return True
    except Exception as e:
        print(f"Error sending command '{command}': {str(e)}")
        return False


def snapshot_regions(bam_file_list, region_list, output_dir, genome='hg38'):
    """
    Generate IGV snapshots in batch.
    """
    sock = None
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory: {os.path.abspath(output_dir)}")

        # Connect to IGV
        sock = connect_igv()

        # Start new session
        send_command(sock, 'new')

        # Set max panel height
        send_command(sock, 'maxPanelHeight 500')

        # Load reference genome
        print(f"Loading genome: {genome}")
        send_command(sock, f'genome {genome}')

        for index, bam_file in enumerate(bam_file_list):
            print(f"\nProcessing file {index + 1}/{len(bam_file_list)}: {bam_file}")

            if not os.path.exists(bam_file):
                print(f"Error: File not found: {bam_file}")
                continue

            if not os.path.exists(bam_file + '.bai'):
                print(f"Warning: Index file not found: {bam_file}.bai")

            if not send_command(sock, f'load {bam_file}'):
                print(f"Failed to load file: {bam_file}")
                continue

            # Key preferences
            send_command(sock, 'preference SAM.SHOW_SMALL_INDELS false')
            send_command(sock, 'preference SAM.SMALL_INDEL_BP_THRESHOLD 50')
            send_command(sock, 'preference SAM.SHOW_INDEL_SIZE true')

            for chrom, start, end, name in region_list:
                try:
                    start, end = start - 1000, end + 1000  # Extend region
                    print(f"\nProcessing region: {chrom}:{start}-{end}")

                    send_command(sock, f'goto {chrom}:{start}-{end}')
                    send_command(sock, 'snapshotSize 800')
                    send_command(sock, 'collapse')

                    snapshot_file = f'{output_dir}/{name}_{index}.png'
                    print(f"Generating snapshot: {snapshot_file}")
                    send_command(sock, f'snapshot {snapshot_file}')

                    if os.path.exists(snapshot_file):
                        print(f"Snapshot generated: {snapshot_file}")
                    else:
                        print(f"Warning: Snapshot not generated: {snapshot_file}")

                except Exception as e:
                    print(f"Error processing region: {str(e)}")
                    continue

            send_command(sock, 'clear')

    except Exception as e:
        print(f"Error occurred: {str(e)}")
        raise
    finally:
        if sock:
            try:
                sock.close()
                print("IGV connection closed")
            except:
                pass


def add_border_and_padding(image, border_size=2, padding=10, border_color='black', bg_color='white'):
    """
    Add border and padding to an image.
    """
    image_with_border = ImageOps.expand(image, border=border_size, fill=border_color)
    final_image = ImageOps.expand(image_with_border,
                                  border=(0, padding // 2, 0, padding // 2),
                                  fill=bg_color)
    return final_image


def combine_images(image_files, output_path, border_size=2, padding=10):
    """
    Combine multiple images vertically into one.
    """
    images = [Image.open(img) for img in image_files]
    images = [add_border_and_padding(img, border_size, padding) for img in images]
    width = max(img.width for img in images)
    total_height = sum(img.height for img in images)

    combined_image = Image.new('RGB', (width, total_height))
    y_offset = 0
    for img in images:
        x_offset = (width - img.width) // 2
        combined_image.paste(img, (x_offset, y_offset))
        y_offset += img.height

    combined_image.save(output_path)
    print(f"Combined image saved: {output_path}")


def process_images(directory):
    """
    Process all PNG files in a directory.
    """
    png_files = [f for f in os.listdir(directory) if f.endswith('.png')]
    groups = collections.defaultdict(list)

    for filename in png_files:
        parts = filename.split('_')
        if len(parts) > 1:
            key = parts[0] + '_' + parts[1]
            full_path = os.path.join(directory, filename)
            groups[key].append(full_path)

    for key, files in groups.items():
        if len(files) > 1:
            files.sort()
            output_path = os.path.join(directory, f'{key}.png')
            try:
                combine_images(files, output_path)
                print(f"Group {key} images combined, total {len(files)} files")
            except Exception as e:
                print(f"Error processing group {key}: {str(e)}")


if __name__ == '__main__':
    output_dir = '/home/huheng/code/PanSV/trio/intersection_eval_hifi/'

    bam_file = [
        '/data1/huheng/HG002/hifi.bam',
        '/data1/huheng/HG003/hg003_hifi.bam',
        '/data1/huheng/HG003/hg004_hifi.bam'
    ]
    regions = []

    print("Reading regions...")
    with open('/home/huheng/code/PanSV/trio/intersect.bed', 'r') as f:
        for line in f:
            chr, start, end, svtype = line.strip().split('\t')
            start, end = int(start), int(end)
            regions.append((chr, start, end, f'{chr}_{start}_{end}'))

    print(f"Total {len(regions)} regions read.")

    snapshot_regions(
        bam_file_list=bam_file,
        region_list=regions,
        output_dir=output_dir,
        genome='hg19'
    )

    process_images(output_dir)
