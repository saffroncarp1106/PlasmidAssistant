def parse_packets(data: bytes) -> list[tuple[int, bytes, bytes]]:
    if len(data) < 5 or data[0] != 9:
        raise ValueError("Invalid first packet header")
    offset = 0
    result = []
    while offset < len(data):
        if offset + 5 > len(data):
            raise ValueError("Truncated header")
        length = int.from_bytes(data[offset+1:offset+5], 'big')
        end = offset + 5 + length
        if end > len(data):
            raise ValueError("Truncated payload")
        packet_type = data[offset]
        payload = data[offset+5:end]
        original_packet = data[offset:end]
        result.append((packet_type, payload, original_packet))
        offset = end
    return result
